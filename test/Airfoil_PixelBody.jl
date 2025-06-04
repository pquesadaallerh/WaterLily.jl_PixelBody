using WaterLily
using Plots

cID = "AirfoilPixelBody"

# TODO: TEMP testing of PixelBody functionality
body = PixelBody(image_to_field("test/resources/airfoil.png")) # n=450, m = 800
img = test_pixel_body(body)
heatmap(img[:,:,1], aspect_ratio=:equal, color=:grays)

# function airfoil(Re=250,U=1,mem=Array)
#     body = PixelBody(image_to_field("test/resources/airfoil.png"))
#     test_pixel_body(body)
#     n = body.n ; m = body.m
#     characteristic_length = m/10 # TODO: Arbitrary atm
#     Simulation((n,m), (U,0), characteristic_length; ν=U*characteristic_length/Re, body, mem) #TODO: Figure out how to estimate chord length from pixel image
#                                                                    # rather than given as an input
# end

# # Initialize the simulation with GPU Array
# # using CUDA
# # sim = airfoil(mem=CuArray);
# sim = airfoil();

# WaterLily.logger(cID) # Log the residual of pressure solver
# #= NOTE: 
# If you want to log residuals during a GPU simulation, it's better to include the following line. 
# Otherwise, Julia will generate excessive debugging messages, which can significantly slow down the simulation. 
# =#
# using Logging; disable_logging(Logging.Debug)

# # Run the simulation
# sim_gif!(sim,duration=10,clims=(-5,5),plotbody=true)

# # Remember to call Plots package (already done in Line 2). This will let WaterLily
# # knows you want to plot sth like residual and will compile the funciton for you.
# # NOTE: Comment out this line if you want to see gif animation!
# # plot_logger("$(cID).log")