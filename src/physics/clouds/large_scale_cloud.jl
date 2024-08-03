abstract type AbstractLargeScaleCloud <: AbstractParameterization end

export AbstractLargeScaleCloud
export NoLargeScaleCloud

struct NoLargeScaleCloud <: AbstractLargeScaleCloud end
NoLargeScaleCloud(::SpectralGrid) = NoLargeScaleCloud()
initialize!(::NoLargeScaleCloud, ::PrimitiveEquation) = nothing
large_scale_cloud!(::ColumnVariables, ::NoLargeScaleCloud, ::PrimitiveEquation) = nothing

# do nothing fall back for primitive dry 
function large_scale_cloud!( 
    column::ColumnVariables,
    model::PrimitiveEquation,
)
    return nothing
end
