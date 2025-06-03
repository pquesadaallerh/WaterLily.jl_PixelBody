"""
    Scale(inS, inE, outS, outE, r)

  - `inS::Float64`: First value of input map range
  - `inE::Float64`: Last value of input map range
  - `outS::Float64`: First value of output map range (i.e. pixel map)
  - `outE::Float64`: Last value of output map range (i.e. pixel map)
  - `r::Float64`: Scaling factor (input to output conversion factor)

Helper struct used together with 'Window' to map between two coordinate systems (i.e. a logical/physical scale 
into a pixel scale). Function "out" converts from logical into pixel coordinates, while "in" converts
from pixel coordinates into logical/physical coordinates.
"""

mutable struct Scale
    inS::Float64
    inE::Float64
    outS::Float64
    outE::Float64
    r::Float64
    
    function Scale(inS::Float64, inE::Float64, outS::Float64, outE::Float64)
        r = (outE - outS) / (inE - inS)
        new(inS, inE, outS, outE, r)
    end
end

# Outer constructor for defualt inS=0.0 and inE=1.0
Scale(outS::Float64, outE::Float64) = Scale(0.0, 1.0, outS, outE)

# Convert from logical coordinates into pixel coordinates
convert_logical_to_pixel(s::Scale, input::Float64) = (input - s.inS) * s.r + s.outS

# Convert from logical coordinates into pixel coordiantes (bounded by domain size)
convert_logical_to_pixel_bounded(s::Scale, input::Float64) = out(s, min(max(input, s.inS), s.inE))

# Convert from pixel coordinates into logical coordinates
convert_pixel_to_logical(s::Scale, output::Float64) = (output - s.outS) / s.r + s.inS


"""
    Window(n0, m0, dn, dm, x0, y0, dx, dy)

Creates a window mapping between a logical grid and a rectangular region of screen pixels.

# Arguments
- `n0::Float64`, `m0::Float64`: Starting logical coordinates (typically cell-centered)
- `dn::Float64`, `dm::Float64`: Number of logical cells in x and y direction
- `x0::Int`, `y0::Int`: Starting pixel coordinates
- `dx::Int`, `dy::Int`: Width and height in pixels

# Fields
- `x::Scale`, `y::Scale`: `Scale` objects mapping logical to pixel coordinates in x and y axes
- `x0::Int`, `y0::Int`: Starting pixel coordinate (same as constructor input)
- `dx::Int`, `dy::Int`: Size in pixels (same as constructor input)

Use this struct to transform between simulation coordinates and display coordinates.
"""

mutable struct Window
    x::Scale
    y::Scale
    x0::Int
    y0::Int
    dx::Int
    dy::Int
    
    function Window(n0::Float64, m0::Float64, dn::Float64, dm::Float64, 
                   x0::Int, y0::Int, dx::Int, dy::Int)
        
        x_scale = Scale(n0, n0 + dn, Float64(x0), Float64(x0 + dx))
        y_scale = Scale(m0, m0 + dm, Float64(y0), Float64(y0 + dy))
        new(x_scale, y_scale, x0, y0, dx, dy)
    end
end

# Default window dimensions # TODO: Not sure how to pass window full size as defaults.
# The 700 was copied from the LilyPad.pde file from LilyPad
const DEFAULT_WIDTH = 700
const DEFAULT_HEIGHT = 700

# Default constructor between 0 and 1 (uses global width/height equivalent)
Window() = Window(0.0, 0.0, 1.0, 1.0, 0, 0, DEFAULT_WIDTH, DEFAULT_HEIGHT)

# Constructor with domain size only - Fixed type conversion
Window(n::Int, m::Int) = Window(1.0, 1.0, Float64(n-2), Float64(m-2), 0, 0, DEFAULT_WIDTH, DEFAULT_HEIGHT)

# Constructor with start + size
Window(n0::Float64, m0::Float64, dn::Float64, dm::Float64) = Window(n0, m0, dn, dm, 0, 0, DEFAULT_WIDTH, DEFAULT_HEIGHT)

#TODO: Constructor with shift in starting indices to allign with cell centers

# Window coordinate conversion methods (round when converting to pixel domain)
ix(w::Window, i::Int) = convert_pixel_to_logical(w.x, Float64(i))
iy(w::Window, i::Int) = convert_pixel_to_logical(w.y, Float64(i))
px(w::Window, i::Float64) = Int(round(convert_logical_to_pixel(w.x, i)))
py(w::Window, i::Float64) = Int(round(convert_logical_to_pixel(w.y, i)))

# Unit conversion methods based on scaling factor
idx(w::Window, i::Int) = Float64(i) / w.x.r
idy(w::Window, i::Int) = Float64(i) / w.y.r
pdx(w::Window, i::Float64) = Int(round(w.x.r * i))
pdy(w::Window, i::Float64) = Int(round(w.y.r * i))

# Method to check if pixel coordinates inside window
inside(w::Window, x::Int, y::Int) = (x ≥ w.x0) && (x ≤ w.x0 + w.dx) && (y ≥ w.y0) && (y ≤ w.y0 + w.dy)