using StaticArrays
using LinearAlgebra
using FileIO
using Images

# Field struct to represent 2D scalar fields
mutable struct Field{T}
    n::Int
    m::Int
    a::Matrix{T}

    # From scalar value (initial fill value)
    function Field{T}(n::Int, m::Int, init_val::T) where T
        new{T}(n, m, fill(init_val, n, m))
    end

    # From matrix
    function Field{T}(a::Matrix{T}) where T
        n, m = size(a)
        new{T}(n, m, a)
    end
end

# Outer constructor to infer type from scalar
Field(n::Int, m::Int, init_val) = Field{typeof(init_val)}(n, m, init_val)
# Outer constructor to infer type from  matrix
Field(a::Matrix{T}) where T = Field{T}(a)

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

# Convert image to fluid-solid Field.
function image_to_field(path::String; solid_val=1.0, fluid_val=0.0, threshold=0.5)
    img = load(path)                          # Load image
    n, m = size(img)                          # Image dimensions
    gray = Gray.(channelview(img))            # Convert to grayscale (if not already)
    
    if ndims(gray) == 3
        gray = gray[1, :, :]                  # Handle color images by picking 1st channel
    end
    
    field_array = zeros(Float64, n, m)
    
    # Fill field: solid where pixel intensity < threshold, fluid otherwise
    for i in 1:n, j in 1:m
        field_array[i, j] = gray[i, j].val < threshold ? solid_val : fluid_val
    end
    
    return Field(n, m, field_array)
end

# PixelBody struct - implements AbstractBody interface
mutable struct PixelBody <: AbstractBody
    pix::Field
    n::Int
    m::Int
    window::Window
    area::Float64
    mass::Float64
    velocity::SVector{2,Float64}  # Changed to SVector for consistency
    
    function PixelBody(pix::Field, velocity::AbstractVector = [0.0, 0.0])
        n, m = pix.n, pix.m
        window = Window(n, m)
        pb = new(pix, n, m, window, 0.0, 0.0, SVector{2,Float64}(velocity))
        get_area!(pb)
        return pb
    end
end

# Constructor with dimensions only
PixelBody(n::Int, m::Int) = PixelBody(Field(n, m), [0.0, 0.0])

# Constructor directly from image (path)
# PixelBody(path::Str) = PixelBody(image_to_field(path), [0.0, 0.0])

# Constructor that creates PixelBody from physical bounds (new functionality)
function PixelBody(pix::Field, x_bounds::Tuple, y_bounds::Tuple, 
                   velocity::AbstractVector = [0.0, 0.0])
    n, m = pix.n, pix.m
    x_start, x_end = x_bounds
    y_start, y_end = y_bounds
    
    # Create window mapping from physical bounds to pixel indices
    window = Window(x_start, y_start, x_end - x_start, y_end - y_start, 
                   1, 1, n, m)
    
    pb = PixelBody.__new__(pix, n, m, window, 0.0, 0.0, SVector{2,Float64}(velocity))
    get_area!(pb)
    return pb
end

# Required AbstractBody interface implementation

"""
    physical_to_grid(body::PixelBody, x)

Convert physical coordinates to grid coordinates using the window mapping.
For backward compatibility, if using default window (grid coordinates), returns x directly.
"""
function physical_to_grid(body::PixelBody, x)
    # Check if using default window (grid coordinates)
    if body.window.x.inS ≈ 1.0 && body.window.y.inS ≈ 1.0 && 
       body.window.x.inE ≈ Float64(body.n-2) && body.window.y.inE ≈ Float64(body.m-2)
        # Default window - assume input is already in grid coordinates
        return [Float64(x[1]), Float64(x[2])]
    else
        # Physical coordinates - convert using window
        return [ix(body.window, px(body.window, x[1])), 
                iy(body.window, py(body.window, x[2]))]
    end
end

"""
    sdf(body::PixelBody, x, t=0)

Signed distance function for PixelBody. Returns the signed distance from point x to the body.
For pixel bodies, this uses the interpolated field values to estimate distance.
"""
function sdf(body::PixelBody, x, t=0.0; kwargs...)
    # Convert physical coordinates to grid coordinates
    grid_coords = physical_to_grid(body, x)
    i, j = grid_coords[1], grid_coords[2]
    
    # Get field value at this point (0 = void, 1 = solid)
    field_val = interp(body.pix, i, j)
    
    # Convert field value to signed distance
    # field_val = 0 means outside (positive distance)
    # field_val = 1 means inside (negative distance)
    # Scale by minimum physical spacing for better accuracy
    min_spacing = min(1.0/body.window.x.r, 1.0/body.window.y.r)
    return (field_val - 0.5) * min_spacing
end

"""
    measure(body::PixelBody, x, t=0; fastd²=Inf)

Returns (d, n, V) where:
- d is the signed distance from x to the body at time t
- n is the normal vector at x
- V is the velocity vector at x

For PixelBody, the normal is computed from the field gradient.
"""
function measure(body::PixelBody, x, t=0.0; fastd² = Inf)
    # Convert physical coordinates to grid coordinates
    grid_coords = physical_to_grid(body, x)
    i, j = grid_coords[1], grid_coords[2]
    
    # Calculate signed distance
    d = sdf(body, x, t)
    
    # Fast approximation if distance is large
    if d^2 > fastd²
        return d, zeros(SVector{2,Float64}), zeros(SVector{2,Float64})
    end
    
    # Calculate normal vector from gradient (using existing dpix function)
    grad = dpix(body, i, j)
    grad_magnitude = norm(grad)
    
    # Convert gradient to physical coordinates
    if grad_magnitude > eps(Float64)
        # Scale gradient by window scaling factors
        grad_physical = [grad[1] * body.window.x.r, grad[2] * body.window.y.r]
        grad_physical_magnitude = norm(grad_physical)
        n = SVector{2,Float64}(grad_physical / grad_physical_magnitude)
    else
        n = zeros(SVector{2,Float64})
    end
    
    # Velocity vector
    V = body.velocity
    
    return d, n, V
end

# Convenience method to set body velocity
function set_velocity!(body::PixelBody, velocity::AbstractVector)
    body.velocity = SVector{2,Float64}(velocity)
end

# Calculate body area
function get_area!(pb::PixelBody)
    pb.area = (pb.n - 2) * (pb.m - 2) - sum(pb.pix)
    pb.mass = pb.area  # default unit density
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

#=TODO: TEMP testing functions for PixelBody =#

# Set rectangular region in field to specified value
function eq!(field::Field, value::Float64, i1::Int, i2::Int, j1::Int, j2::Int)
    for i in i1:i2, j in j1:j2
        if 1 <= i <= field.n && 1 <= j <= field.m
            field.a[i, j] = value
        end
    end
end

# Example usage and setup functions
function create_example_body()
    n = Int(2^7)  # 128 x 128 grid points
    field = Field(n, n, 1.0)  # Initialize with 1.0 (fluid)
    
    # Add void regions (set to 0.0)
    # eq!(field, 0.0, 40, 50, 30, 40)   # void square
    eq!(field, 0.0, 20, 70, 60, 65)   # void rectangle  
    # eq!(field, 0.0, 40, 50, 80, 90)   # void square

    return PixelBody(field)
end

# Display function (returns RGB array for visualization)
function display_data(pb::PixelBody)
    img_data = zeros(Float64, pb.window.dy, pb.window.dx, 3)  # RGB array
    
    for i in 1:pb.window.dx
        x = ix(pb.window, i + pb.window.x0)
        for j in 1:pb.window.dy
            y = iy(pb.window, j + pb.window.y0)
            f = interp(pb.pix, x, y)
            
            # Create grayscale visualization where f=1 is white, f=0 is black
            gray_val = f
            img_data[j, i, 1] = gray_val  # R
            img_data[j, i, 2] = gray_val  # G  
            img_data[j, i, 3] = gray_val  # B
        end
    end
    
    return img_data
end

# Helper function to visualize the body (for testing)
function test_pixel_body(body::PixelBody)
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
    
    return img_data
end

# Usage example:
# body = create_example_body()
# img = test_pixel_body(body)
# using Plots
# heatmap(img[:,:,1], aspect_ratio=:equal, color=:grays)
