
using PyPlot
include("../util/misc.jl")
"""
    animated_gif(stimulus...; filename="default.gif", verbose=0, kwargs...)

Creates an animated gif of the given stimuli (can animate more than one at a
time for side-by-side comparison). Super hack-y: Uses `PyPlot`'s `imshow`
function for each frame of the stimulus, writes them to separate `.gif` files,
then uses a system call to `convert` to create the animated gif. Will delete the
individual frames after it's done, and will do its best to be graceful about
this (i.e. if it's interrupted, that sort of thing).

Note that it will remove any filename extensions from `filename` using
`remove_extension` and then append ".gif". In particular, if there are any `'.'`
characters in `filename`, make sure it ends with `.gif`

Returns `success, status` where `success` is a boolean indicated success or
failure of the process to produce an animated gif, and `status` is an array of
strings representing status messages produced throughout (e.g. warnings about
which frames got skipped, that sort of thing).

"""
function animated_gif(stimulus...; filename="default.gif", verbose=0,
    kwargs...)

    dkwargs = Dict(kwargs)
    floc = dirname(abspath(filename))
    fname_root = remove_extension(basename(filename))
    temp_dir = joinpath(floc, fname_root, "temp")
    gif_filename = joinpath(floc, fname_root * ".gif")

    image_arrays = [frame_image(S, frame_range) for S in stimulus]

    status = String[]

    start_frame = pop!(dkwargs, :start_frame, 1)
    end_frame = pop!(dkwargs, :end_frame, 10)
    frame_range = pop!(dkwargs, :frame_range, start_frame:end_frame)

    cmap = pop!(dkwargs, :cmap, "gray")
    aspect = pop!(dkwargs, :aspect, "equal")
    cbar = pop!(dkwargs, :colorbar, true)
    normalize = pop!(dkwargs, :normalize, false)
    if normalize
        vmin = fill(minimum(map(minimum,image_arrays)), size(image_arrays))
        vmax = fill(maximum(map(maximum,image_arrays)), size(image_arrays))
    else
        vmin = map(minimum, image_arrays)
        vmax = map(maximum, image_arrays)
    end

    layout_x = ceil(Int, sqrt(length(stimulus)))
    layout_y = ceil(Int, sqrt(length(stimulus)))

    subplot_titles = pop!(dkwargs, :subplot_titles, ["Stimulus $i" for i in 1:length(stimulus)])

    try
        mkpath(temp_dir)

        for frame_idx in 1:length(frame_range)
            frame_filename = "$fname_root-" * @sprintf("%04d",idx) * ".gif"

            fig = figure("Frame $frame_idx")
            for stim_idx in 1:length(image_arrays)
                subplot("$layout_x$layout_y$stim_idx")
                imshow(image_arrays[stim_idx][:,:,frame_idx], cmap=cmap, aspect=aspect, vmin=vmin[stim_idx], vmax=vmax[stim_idx])
                if cbar
                    colorbar()
                end
                title(subplot_titles[stim_idx])
                # TODO finish this whole thing up, still a lot to be done but it can wait
            end
            savefig(joinpath(temp_dir, frame_filename))
            close(fig)
        end
        # for (idx, frame_number) in enumerate(frame_range)
        #     # if, for whatever reason, frame_range is passed with some weird order,
        #     # I want the resulting gif to respect that order.
        #     frame_filename = "$fname_root-" * @sprintf("%04d",idx) * ".gif"
        #
        #     fig = figure()
        #     for (idx, S) in enumerate(stimulus)
        #         subplot("$layout_x$layout_y$idx")
        #         imshow(frame_image(S, frame_number), cmap=cmap, aspect=aspect)
        #         if colorbar
        #             colorbar()
        #         end
        #     end
        #     close(fig)
        # end
        success = true
    catch except
        success = false
    finally
        close("all")
        rm(joinpath(floc, fname_root), force=true, recursive=true)
    end
    return success, status
end
