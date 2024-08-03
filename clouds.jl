### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 4e079cea-492f-11ef-2b79-2f7a0c3eb938
begin
	import Pkg
	Pkg.activate(Base.current_project())
	# Pkg.build()
	Pkg.instantiate()
	# Pkg.test()
	using Revise
	using CairoMakie
	using SpeedyWeather
	import CUDA
	# WARNING TO SELF: "Revise" cannot deal with modified kwdef in structs. If you modify one of those you probably need to restart the notebook.
end

# ╔═╡ fead1827-28bc-4298-b8ac-6344dec5e719
spectral_grid = SpectralGrid(trunc=50, nlev=8)

# ╔═╡ 1b02f213-cbbe-412b-a30e-0cd8026d1f21
# ╠═╡ disabled = true
#=╠═╡
condensation = CloudsWithImplicitCondensation(spectral_grid, ls_cloud_scheme=HomogeneousCloud(spectral_grid, cloud_cover_scheme=MaximumRandomOverlap(spectral_grid)))
  ╠═╡ =#

# ╔═╡ c67c36c3-757e-42b5-b911-332e28a6ab47
condensation = CloudsWithImplicitCondensation(spectral_grid, ls_cloud_scheme=SmithCloud(spectral_grid, cloud_cover_scheme=MaximumRandomOverlap(spectral_grid)))

# ╔═╡ 9c58d48f-47cd-49ba-a041-20cc972a6302
# ╠═╡ disabled = true
#=╠═╡
condensation = ImplicitCondensation(spectral_grid)
  ╠═╡ =#

# ╔═╡ 16c9c4b2-2be7-4335-9882-d4e76bc5b57e
model = PrimitiveWetModel(; spectral_grid, large_scale_condensation=condensation, shortwave_radiation=NoShortwave())

# ╔═╡ 74db0f09-d22a-440b-a0ec-767fb5eb1244
simulation = SpeedyWeather.initialize!(model)

# ╔═╡ 409579a5-8a1d-499e-a96b-606845d38d9b
run!(simulation, period=Day(60))

# ╔═╡ 37b10060-0a70-43a8-9c3e-c8079cbc7fe1
# ╠═╡ disabled = true
#=╠═╡
heatmap(simulation.diagnostic_variables.surface.cloud_top)
  ╠═╡ =#

# ╔═╡ e691d99d-ce62-40de-aa9d-75286f90bd9d
heatmap(simulation.diagnostic_variables.surface.cloud_cover)

# ╔═╡ 19b58cfb-a589-4dfa-8ae4-b5bcfdf9bc80
heatmap(simulation.diagnostic_variables.surface.cloud_top)

# ╔═╡ 7bdbdad5-e486-449d-9537-8193b6163280
heatmap(simulation.diagnostic_variables.surface.precip_large_scale)

# ╔═╡ cb04f1df-bd40-4f37-9405-eea171e122e3
heatmap(simulation.diagnostic_variables.surface.precip_convection)

# ╔═╡ 13e5a1f3-5f87-4ae4-af6c-2cc2bf24a32a
heatmap(simulation.diagnostic_variables.surface.precip_convection + 
simulation.diagnostic_variables.surface.precip_large_scale)

# ╔═╡ 9fafa523-0f19-48cd-9754-c5dba168dd32
CUDA.functional()

# ╔═╡ 61711c12-a759-4c98-8b4b-8de59dcdcad1
spectral_grid2 = SpectralGrid(trunc=85, nlev=10)

# ╔═╡ c4097121-0db9-414b-8521-39c8f831a1db
condensation2 = CloudsWithImplicitCondensation(spectral_grid2, ls_cloud_scheme=SlingoCloud(spectral_grid2, cloud_cover_scheme=MaximumRandomOverlap(spectral_grid2)));

# ╔═╡ e998505a-da55-40b9-9d08-0051dea7787b
model2 = PrimitiveWetModel(; spectral_grid=spectral_grid2, large_scale_condensation=condensation2);

# ╔═╡ 2d94e9b0-42a1-4151-bf6b-9038c857812e
simulation2 = SpeedyWeather.initialize!(model2);

# ╔═╡ f22a3a68-9f97-48a8-9073-e328bc9ca323
run!(simulation2, period=Day(120))

# ╔═╡ f925e22a-577c-4006-8310-5e61a2542310
heatmap(simulation2.diagnostic_variables.surface.cloud_cover)

# ╔═╡ 344642ac-9dda-4cd0-8c57-75e22bdf6a44
heatmap(simulation2.diagnostic_variables.surface.cloud_top)

# ╔═╡ f7c6363e-ad56-4c8d-83a3-b37754c7c7cf
heatmap(simulation2.diagnostic_variables.surface.precip_large_scale)

# ╔═╡ 057508c7-ec7b-404f-bf61-d5d324c44cf2
heatmap(simulation2.diagnostic_variables.surface.precip_convection)

# ╔═╡ 8a1a2ad5-9e56-496e-8969-61412a71cd88
heatmap(simulation2.diagnostic_variables.surface.precip_convection + 
simulation2.diagnostic_variables.surface.precip_large_scale)

# ╔═╡ Cell order:
# ╠═4e079cea-492f-11ef-2b79-2f7a0c3eb938
# ╠═fead1827-28bc-4298-b8ac-6344dec5e719
# ╠═1b02f213-cbbe-412b-a30e-0cd8026d1f21
# ╠═c67c36c3-757e-42b5-b911-332e28a6ab47
# ╠═9c58d48f-47cd-49ba-a041-20cc972a6302
# ╠═16c9c4b2-2be7-4335-9882-d4e76bc5b57e
# ╠═74db0f09-d22a-440b-a0ec-767fb5eb1244
# ╠═409579a5-8a1d-499e-a96b-606845d38d9b
# ╠═37b10060-0a70-43a8-9c3e-c8079cbc7fe1
# ╠═e691d99d-ce62-40de-aa9d-75286f90bd9d
# ╠═19b58cfb-a589-4dfa-8ae4-b5bcfdf9bc80
# ╠═7bdbdad5-e486-449d-9537-8193b6163280
# ╠═cb04f1df-bd40-4f37-9405-eea171e122e3
# ╠═13e5a1f3-5f87-4ae4-af6c-2cc2bf24a32a
# ╠═9fafa523-0f19-48cd-9754-c5dba168dd32
# ╠═61711c12-a759-4c98-8b4b-8de59dcdcad1
# ╠═c4097121-0db9-414b-8521-39c8f831a1db
# ╠═e998505a-da55-40b9-9d08-0051dea7787b
# ╠═2d94e9b0-42a1-4151-bf6b-9038c857812e
# ╠═f22a3a68-9f97-48a8-9073-e328bc9ca323
# ╠═f925e22a-577c-4006-8310-5e61a2542310
# ╠═344642ac-9dda-4cd0-8c57-75e22bdf6a44
# ╠═f7c6363e-ad56-4c8d-83a3-b37754c7c7cf
# ╠═057508c7-ec7b-404f-bf61-d5d324c44cf2
# ╠═8a1a2ad5-9e56-496e-8969-61412a71cd88
