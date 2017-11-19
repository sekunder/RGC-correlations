
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
strings and `Exception`s representing status messages produced throughout (e.g.
warnings about which frames got skipped, that sort of thing).

"""
function animated_gif(stimulus...; filename="default.gif", verbose=0,
    kwargs...)

    ########################################
    #### Handle filename stuff
    ########################################
    floc = dirname(abspath(filename))
    fname_root = remove_extension(basename(filename))
    temp_dir = joinpath(floc, fname_root, "temp")
    gif_filename = joinpath(floc, fname_root * ".gif")

    dkwargs = Dict(kwargs)

    ########################################
    #### Get animation properties
    ########################################
    start_frame = pop!(dkwargs, :start_frame, 1)
    end_frame = pop!(dkwargs, :end_frame, 10)
    frame_range = pop!(dkwargs, :frame_range, start_frame:end_frame)
    loop = pop!(dkwargs, :loop, 0)
    # Framerate: in order, check for :fps, :frame_time_s
    fps = pop!(dkwargs, :fps, 60)
    frame_time_s = pop!(dkwargs, :frame_time_s, round(Int,1/fps))
    frame_time_hundredths = round(Int,100frame_time_s)

    ########################################
    #### Get pyplot options
    ########################################
    cmap = pop!(dkwargs, :cmap, "gray")
    aspect = pop!(dkwargs, :aspect, "equal")
    cbar = pop!(dkwargs, :colorbar, true)
    tight = pop!(dkwargs, :tight_layout, true)
    normalize = pop!(dkwargs, :normalize, false)
    if normalize
        # normalize across stimuli as well as across frames
        vmin = fill(minimum(map(minimum,image_arrays)), size(image_arrays))
        vmax = fill(maximum(map(maximum,image_arrays)), size(image_arrays))
    else
        # normalize frames so that a given stimulus uses the same range for each frame.
        vmin = map(minimum, image_arrays)
        vmax = map(maximum, image_arrays)
    end
    layout_x = ceil(Int, sqrt(length(stimulus)))
    layout_y = ceil(Int, sqrt(length(stimulus)))
    titles = pop!(dkwargs, :titles, ["Stimulus $i" for i in 1:length(stimulus)])

    ########################################
    #### State setting things up
    ########################################
    # TODO this is dangerous when frame_range is large!
    image_arrays = [frame_image(S, frame_range) for S in stimulus]

    status = []

    try
        mkpath(temp_dir)

        # for frame_idx in 1:length(frame_range)
        for (frame_idx, frame_number) in enumerate(frame_range)
            frame_filename = "$fname_root-" * @sprintf("%04d",idx) * ".gif"

            fig = figure("Frame $frame_idx", tight_layout=true)
            # for stim_idx in 1:length(image_arrays)
            for (stim_idx,(S,img_arr)) in enumerate(zip(stimulus,image_arrays))
                subplot("$layout_x$layout_y$stim_idx")
                imshow(img_arr[:,:,frame_idx], cmap=cmap, aspect=aspect, vim=vmin[stim_idx], vmax=vmax[stim_idx])
                if cbar
                    colorbar()
                end
                title(titles[stim_idx])
                # TODO could add more options for xticks, yticks, that sort of thing later

                # for now, clear ticks
                xticks()
                yticks()

                # for now, xlabel is "t = X s"
                xlabel("t = " * @sprintf("%0.3f", index_to_time(frame_idx)) * " s")
            end
            savefig(joinpath(temp_dir, frame_filename))
            close(fig)
        end

        # Now, let's do a system call to use convert
        convert_cmd = `convert -delay $frame_time_hundredths -loop $loop $(joinpath(temp_dir,fname_root))-*.gif $gif_filename`
        run(convert_cmd)
        success = true
    catch except
        success = false
        push!(status, except)
    finally
        close("all")
        rm(joinpath(floc, fname_root), force=true, recursive=true)
    end
    return success, status
end
