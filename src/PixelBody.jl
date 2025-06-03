using LinearAlgebra

# Field struct to represent 2D scalar fields
mutable struct Field
    n::Int
    m::Int
    a::Matrix{Float64}
    
    function Field(n::Int, m::Int, init_val::Float64 = 0.0)
        new(n, m, fill(init_val, n, m))
    end
end

# Set rectangular region in field to specified value
function eq!(field::Field, value::Float64, i1::Int, i2::Int, j1::Int, j2::Int)
    for i in i1:i2, j in j1:j2
        if 1 <= i <= field.n && 1 <= j <= field.m
            field.a[i, j] = value
        end
    end
end

# Sum all values in the field
function Base.sum(field::Field)
    return sum(field.a)
end

# Bilinear interpolation for field values
function interp(field::Field, x::Float64, y::Float64)
    # Clamp coordinates to field bounds
    x = clamp(x, 1.0, Float64(field.n))
    y = clamp(y, 1.0, Float64(field.m))
    
    i = Int(floor(x))
    j = Int(floor(y))
    
    # Handle boundary cases
    i = clamp(i, 1, field.n-1)
    j = clamp(j, 1, field.m-1)
    
    # Fractional parts
    fx = x - i
    fy = y - j
    
    # Bilinear interpolation
    c00 = field.a[i, j]
    c10 = field.a[i+1, j]
    c01 = field.a[i, j+1]
    c11 = field.a[i+1, j+1]
    
    c0 = c00 * (1 - fx) + c10 * fx
    c1 = c01 * (1 - fx) + c11 * fx
    
    return c0 * (1 - fy) + c1 * fy
end

# PixelBody struct - implements AbstractBody interface
mutable struct PixelBody <: AbstractBody
    pix::Field
    n::Int
    m::Int
    window::Window
    area::Float64
    mass::Float64
    velocity::Vector{Float64}  # Body velocity for measure function
    
    function PixelBody(pix::Field, velocity::Vector{Float64} = [0.0, 0.0])
        n, m = pix.n, pix.m
        window = Window(n, m)
        pb = new(pix, n, m, window, 0.0, 0.0, velocity)
        get_area!(pb)
        return pb
    end
end

# Constructor with dimensions only
PixelBody(n::Int, m::Int) = PixelBody(Field(n, m), [0.0, 0.0])

# Required AbstractBody interface implementation

"""
    sdf(body::PixelBody, x, t=0)

Signed distance function for PixelBody. Returns the signed distance from point x to the body.
For pixel bodies, this uses the interpolated field values and gradients to estimate distance.
"""
function sdf(body::PixelBody, x::AbstractVector, t::Real = 0.0; kwargs...)
    # Convert physical coordinates to grid coordinates
    i = Float64(x[1])
    j = Float64(x[2])
    
    # Get field value at this point (0 = void, 1 = solid)
    field_val = interp(body.pix, i, j)
    
    # Convert field value to signed distance
    # field_val = 0 means outside (positive distance)
    # field_val = 1 means inside (negative distance)
    # Use a simple linear mapping for now
    return field_val - 0.5  # Maps [0,1] to [-0.5, 0.5]
end

"""
    measure(body::PixelBody, x, t=0, fastd²=Inf)

Returns (d, n, V) where:
- d is the signed distance from x to the body at time t
- n is the normal vector at x
- V is the velocity vector at x

For PixelBody, the normal is computed from the field gradient.
"""
function measure(body::PixelBody, x; fastd² = Inf)
    # Convert physical coordinates to grid coordinates
    i = Float64(x[1])
    j = Float64(x[2])
    
    # Calculate signed distance
    d = sdf(body, x, t)
    
    # Fast approximation if distance is large
    if d^2 > fastd²
        return d, zero(x), zero(x)
    end
    
    # Calculate normal vector from gradient
    grad = dpix(body, i, j)
    n = length(grad) > 0 ? normalize(grad) : zero(x)
    
    # Velocity vector (constant for static bodies, could be time-dependent)
    V = copy(body.velocity)
    
    return d, n, V
end

# Convenience method to set body velocity
function set_velocity!(body::PixelBody, velocity::Vector{Float64})
    body.velocity = velocity
end

# Calculate body area
function get_area!(pb::PixelBody)
    pb.area = (pb.n - 2) * (pb.m - 2) - sum(pb.pix)
    pb.mass = pb.area  # default unit density
end

# Display function (returns RGB array for visualization)
function display_data(pb::PixelBody)
    img_data = zeros(Float64, pb.window.dy, pb.window.dx, 3)  # RGB array
    
    for i in 1:pb.window.dx
        x = ix(pb.window, i + pb.window.x0)
        for j in 1:pb.window.dy
            y = iy(pb.window, j + pb.window.y0)
            f = interp(pb.pix, x, y)
            
            # Create grayscale visualization where f=0 is white, f=1 is black
            gray_val = 1.0 - f
            img_data[j, i, 1] = gray_val  # R
            img_data[j, i, 2] = gray_val  # G  
            img_data[j, i, 3] = gray_val  # B
        end
    end
    
    return img_data
end

# Gradient of pixel field at point (i,j)
function dpix(pb::PixelBody, i::Float64, j::Float64)
    dx = interp(pb.pix, i + 0.5, j) - interp(pb.pix, i - 0.5, j)
    dy = interp(pb.pix, i, j + 0.5) - interp(pb.pix, i, j - 0.5)
    return [dx, dy]
end

# Distance function with gradient magnitude
function del(pb::PixelBody, i::Float64, j::Float64, eps::Float64)
    grad = dpix(pb, i, j)
    return [interp(pb.pix, i, j), norm(grad) * eps]
end

# Wall normal vector
function wall_normal(pb::PixelBody, i::Float64, j::Float64)
    n = dpix(pb, i, j)
    mag = norm(n)
    if mag > 1e-8
        n = n / mag
    end
    return n
end

# Compute pressure force on PixelBody
function press_force(pb::PixelBody, p::Field)
    pv = [0.0, 0.0]
    
    for i in 2:(pb.n-1)
        for j in 2:(pb.m-1)
            grad = dpix(pb, Float64(i), Float64(j))
            pv += grad * p.a[i, j]
        end
    end
    
    return pv
end

# Example usage and setup functions
function create_example_body()
    n = Int(2^7)  # 128 grid points
    field = Field(n, n, 1.0)  # Initialize with 1.0 (solid)
    
    # Add void regions (set to 0.0)
    eq!(field, 0.0, 40, 50, 30, 40)   # void square
    eq!(field, 0.0, 50, 60, 50, 70)   # void rectangle  
    eq!(field, 0.0, 40, 50, 80, 90)   # void square
    
    return PixelBody(field)
end

# Helper function to visualize the body (for testing)
function test_pixel_body()
    body = create_example_body()
    img_data = display_data(body)
    
    println("Created PixelBody with dimensions: $(body.n) x $(body.m)")
    println("Body area: $(body.area)")
    println("Body mass: $(body.mass)")
    println("Image data shape: $(size(img_data))")
    
    # Test AbstractBody interface
    test_point = [45.0, 35.0]
    d = sdf(body, test_point)
    d_measure, n, V = measure(body, test_point)
    
    println("SDF at $test_point: $d")
    println("Measure at $test_point: d=$d_measure, n=$n, V=$V")
    
    # Test gradient calculation
    grad = dpix(body, 45.0, 35.0)
    println("Gradient at (45, 35): $grad")
    
    # Test wall normal
    normal = wall_normal(body, 45.0, 35.0)
    println("Wall normal at (45, 35): $normal")
    
    return body, img_data
end


# Usage example:
# body, img = test_pixel_body()
# using Plots
# heatmap(img[:,:,1], aspect_ratio=:equal, color=:grays)
