abstract type CloudCover <: AbstractParameterization end
# --- Cloud cover schemes ---
export CloudCover
export RandomOverlap
export MaximumRandomOverlap
export MaximumOverlap

struct RandomOverlap <: CloudCover end
RandomOverlap(::SpectralGrid) = RandomOverlap()
initialize!(::RandomOverlap, ::PrimitiveEquation) = nothing
function cloud_cover!(column::ColumnVariables, scheme::RandomOverlap)
    column.cloud_cover = 1
    @inbounds for k in eachindex(column)
        column.cloud_cover *= (1 - column.bulk_cloud_fraction[k])
    end
    column.cloud_cover = 1 - column.cloud_cover
end

struct MaximumRandomOverlap <: CloudCover end
MaximumRandomOverlap(::SpectralGrid) = MaximumRandomOverlap()
initialize!(::MaximumRandomOverlap, ::PrimitiveEquation) = nothing
function cloud_cover!(column::ColumnVariables, scheme::MaximumRandomOverlap)
    column.cloud_cover = 1
    @inbounds for k in eachindex(column)
        # continue where k-1 is out of bounds
        if k == 1
            continue
        end
        # if the cloud fraction is 1 in level k-1, break and return 0
        if isapprox(column.bulk_cloud_fraction[k-1], 1)
            column.cloud_cover = 0
            break
        end
        column.cloud_cover *= (1-max(column.bulk_cloud_fraction[k], column.bulk_cloud_fraction[k-1])/(1-column.bulk_cloud_fraction[k-1]))
    end
    column.cloud_cover = 1-column.cloud_cover
    # bad division can cause errors here.
    # bounds check: if cloud cover mra is above 1, set it to 1; if below 0, set it to 0
    column.cloud_cover = min(1, max(0, column.cloud_cover))
end

struct MaximumOverlap <: CloudCover end
MaximumOverlap(::SpectralGrid) = MaximumOverlap()
initialize!(::MaximumOverlap, ::PrimitiveEquation) = nothing
function cloud_cover!(column::ColumnVariables, scheme::MaximumOverlap)
    column.cloud_cover = 0
    @inbounds for k in eachindex(column)
        if column.bulk_cloud_fraction[k] > column.cloud_cover
            column.cloud_cover = column.bulk_cloud_fraction[k]
        end
    end
end
    