
"""
    GrayScaleStimulus(pixel_vals, N, px, d, mm_per_px, frame_length_s, frame_rate, onset, zerotonegative, metadata)

Object that represents a stimulus that is displayed on screen to the retina.
Values in `pixel_vals` should be `Float`s (32 or 64) between 0 and 1, inclusive, or `UInt8`s between 0 and 255, inclusive; `pixel_vals` must be a
matrix with second dimension equal to the number of frames.

The values `N`, `px`, and `d` determine how the image is displayed on screen.
Each should be given in "real world" coordinates, where the first coordinate is
the horizontal one. Then, the fact that the first coordinate of a matrix is the
"vertical" coordinate for the purposes of things like `matshow` from `PyPlot` is
handled within the methods associated to this type, by reversing the tuples as
appropriate.

`zerotonegative` is a `Bool` that indicates if the pixel value range should be
stretched to [-1,1] when performing computations.

"""
type GrayScaleStimulus{T<:AbstractArray} <: AbstractStimulus
    # the actual values for the screen
    pixel_vals::T
    # the dimensions on screen
    N::Vector{Int} # resolution of the image (number of distinct patches)
    px::Vector{Int} # the number of pixels in the full image
    d::Vector{Int} # size of each patch
    mm_per_px::Float64
    # timing. frame_rate = 1/frame_length_s
    frame_length_s::Float64
    frame_rate::Float64
    onset::Float64
    # metadata
    zerotonegative::Bool # whether a 0 in pixel_values should be changed to -1
    metadata::Dict
end

"""
    GrayScaleStimulus(pixel_values, existing_GrayScaleStimulus; kwargs...)

Copies the measurements from `existing_GrayScaleStimulus` into a new `GrayScaleStimulus`
object, using the new pixel values. Overwrites old metadata with kwargs.
"""
function GrayScaleStimulus(values::Union{BitMatrix,Matrix{Float64}}, S::GrayScaleStimulus;
    onset=S.onset, zerotonegative=S.zerotonegative, kwargs...)

    return GrayScaleStimulus(values, S.N, S.px, S.d, S.mm_per_px, S.frame_length_s, S.frame_rate, onset, zerotonegative, Dict(kwargs))
end

function show(io::IO, S::GrayScaleStimulus)
    println(io, "Grayscale stimulus")
    println(io, "Duration: $(frame_time(S) * n_frames(S)) s ($(n_frames(S)) frames)")
    println(io, "Frame size (w,h): $(frame_size(S)) pixels, $(S.mm_per_px .* frame_size(S)) mm")
    println(io, "Resolution (w,h): $(S.N)")
    println(io, "Frame rate: $(S.frame_rate) Hz ($(frame_time(S)) s/frame)")
    show_metadata(io, S)
end

frame_size(S::GrayScaleStimulus) = S.px
frame_time(S::GrayScaleStimulus) = S.frame_length_s
n_frames(S::GrayScaleStimulus) = size(S.pixel_vals, 2)
resolution(S::GrayScaleStimulus) = S.N

"""
    time_to_index(S, t, relative_time=false)

The index of the frame that is on screen at time `t`. If `relative_time=true`,
interprets `t` as seconds since `S.onset`, otherwise interprets as absolute
time.

"""
time_to_index(S::GrayScaleStimulus, t::Float64, relative_time::Bool=false) = 1 + floor(Int, (relative_time ? t : t - S.onset)/frame_time(S))

"""
    index_to_time(S, idx, relative_time=false)

The time at which frame number `idx` is displayed on screen, with the first
frame appearing at time `S.onset` (or `0.0` if `relative_time=true`)

"""
index_to_time(S::GrayScaleStimulus, idx::Integer, relative_time::Bool=false) = (idx - 1) * frame_time(S) + (relative_time ? 0.0 : S.onset)

"""
    matrix_form(S)::Matrix{Float64}

Returns a minimal representation of `S` for performing computations.
"""
matrix_form(S::GrayScaleStimulus, frames=1:size(S.pixel_vals,2)) = _pixel_values_to_float(S.pixel_vals[:,frames], S.zerotonegative)

_pixel_values_to_float(v::BitArray, negative::Bool) = negative ? (-1.0) .^ (.!v) : Matrix{Float64}(v)
_pixel_values_to_float(v::Array{UInt8}, negative::Bool) = _pixel_values_to_float(v / 255.0, negative)
_pixel_values_to_float(v, negative::Bool) = negative ? 2.0 * v .- 1.0 : Matrix{Float64}(v)

"""
    frame_image(S, indexes)
    frame_image(S, time, relative=false)

Returns the frame displayed on screen at the given time. If `relative` is true,
the time is time since stimulus onset in seconds. The first form can accept any
valid way of indexing an array (that is, it will pull up
`S.pixel_vals[:,indexes]`); if `indexes` is a single number, will squeeze the
last dimension (i.e. return a 2d array); otherwise returns a 3d array where the
last index is the frame number.

"""
function frame_image(S::GrayScaleStimulus, indexes)
    if isa(indexes, Number)
        indexes = [indexes]
    end
    img = zeros(reverse(S.px)..., length(indexes))
    for (img_idx, stim_idx) in enumerate(indexes)
        fm = reshape(S.pixel_vals[:,stim_idx], reverse(S.N)...)
        frame = _pixel_values_to_float(fm, S.zerotonegative)
        img[:,:,img_idx] = kron(frame, ones(reverse(S.d)...))
    end
    # fm = reshape(S.pixel_vals[:,idx], reverse(S.N))
    # frame = _pixel_values_to_float(fm, S.zerotonegative)
    # return kron(frame, ones(reverse(S.d)...))
    if length(indexes) == 1
        img = squeeze(img, 3)
    end
    return img
end
frame_image(S::GrayScaleStimulus, t::Float64, relative_time::Bool=false) = frame_image(S, ceil(Int, (relative_time ? t - S.onset : t)/frame_time(S)))

"""
    compute_STRFs(spike_hist, S; kwargs...)

Returns an array of `GrayScaleStimulus` objects, one for each neuron.
"""
function compute_STRFs(spike_hist::Matrix{Float64}, S::GrayScaleStimulus; kwargs...)
    RFs = _compute_STRFs(spike_hist, S; kwargs...)
    # let's pop kwargs that are related to computing the RFs
    dkwargs = Dict{Any,Any}(kwargs)
    window_length_s = pop!(dkwargs, :window_length_s, 0.5)
    # dkwargs[:autocomment] = "STRF computed with compute_STRFs"
    return [GrayScaleStimulus(RFs[:,i,:], S; onset=-window_length_s, zerotonegative=false, autocomment="STRF computed with compute_STRFs", dkwargs...) for i = 1:size(RFs,2)]
end
