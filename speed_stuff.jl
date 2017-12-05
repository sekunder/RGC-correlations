println("-" ^ 80)
println("""
This script demonstrates a funny thing that can happen when processing a list of
numbers. Namely, if the operation you carry out depends on whether the number is
high or low, you can achieve a significant speedup by sorting the list first.

The function we use is very simple: Given a list, it loops through the list and
counts the number of `high` and `low` entries (determined by a threshold which
defaults to 0.5, since we'll be using `rand(n)` as our list -- these values lie
in the half-open interval [0,1)). Here's the definition of the function:

    function count_high_low(L, theta=0.5)
        high,low = 0,0
        for x in L
            if x > theta
                high += 1
            else
                low +=1
            end
        end
        return high, low
    end
""")

function count_high_low(L, theta=0.5)
    high,low = 0,0
    for x in L
        if x > theta
            high += 1
        else
            low +=1
        end
    end
    return high, low
end


# let's for julia to compile this function ahead of time so it doesn't mess with
# our timing
L = rand(5)
L_high, L_low = count_high_low(L)
println("For example, given the list L = $(join(L,", ")),")
println("There are $L_high \"high\" entries")
println("      and $L_low \"low\" entries.")

# Pro tip: This block of code is "initialization" stuff. It doesn't "do"
# anything other than define some values. This script went through a series of
# changes: First I had just one loop, then I added stuff to that one loop, then
# I added a second loop later. To make my life easy, I used this initialization
# block to avoid having to hardcode values in my loops.
n_trials = 10^4
n_rand = 10^5; rand_line = @__LINE__
unsorted_time = 0.0
sorted_time = 0.0
time_to_sort = 0.0
println("Beginning timing. Running $n_trials trials on lists of length $n_rand")
# The loop below is the part of my script that "does stuff". Notice there are
# *no* hardcoded values in here! Sometimes avoiding a hardcoded value is
# difficult/not worth the time to fix. But, whenever possible, try to define a
# variable that stores, say, the number of trials you're going to run or
# whatever else. This way if you add more things later (like I added the loop
# below), or use that value in multiple places in the code, or want to add a
# feature that allows you to change the value interactively, it's easier to do.
for trial = 1:n_trials
    R = rand(n_rand)
    tic()
    count_high_low(R)
    unsorted_time += toq()

    tic()
    sort!(R)
    time_to_sort += toq()

    tic()
    count_high_low(R)
    sorted_time += toq()
end
println("Approximate total run time: $(unsorted_time + time_to_sort + sorted_time) s")

println()
println("Average time to run count_high_low:                " * (@sprintf "%0.10f s" (unsorted_time/n_trials)) * " (a)")
println("Average time to sort the list:                     " * (@sprintf "%0.10f s" (time_to_sort/n_trials)))
println("Average time to run count_high_low after sorting:  " * (@sprintf "%0.10f s" (sorted_time/n_trials)) * " (b)")
println("Average time to sort + run count_high_low:         " * (@sprintf "%0.10f s" ((sorted_time + time_to_sort) / n_trials)) )
println()
println("Ratio of averages (unsorted / sorted): (a) / (b) =")
println("                                                   *** $(unsorted_time / sorted_time) ***")
println()

sort_speedup = unsorted_time / sorted_time
s = round(Int, sort_speedup)
println("""
Notice that, on average, the process runs nearly $s times faster on sorted data!
This has to do with "pipelining" and "branchpoints" and some other technical
nonsense. In this particular example, the time to sort heavily outweighs the
speedup, but that's only because the operation we're doing is so simple
(incrementing a counter).

The deep technical stuff isn't really important. What is important is to stop
and think about what this function is actually doing. If we just want to know
how many entries are above and below a certain threshold, it would be much
faster to look at a sorted list and just find the first entry larger than the
threshold. So let's write a function that does that and see if we get a
noticeable speed up.

We define the following function:

    function count_high_low_by_sort(L, theta=0.5)
        L = sort(L)
        low = searchsortedlast(L, theta)
        return length(L) - low, low
    end
""")

function count_high_low_by_sort(L, theta=0.5)
    L = sort(L)
    low = searchsortedlast(L, theta)
    return length(L) - low, low
end

L_high, L_low = count_high_low_by_sort(L)
println("First, let's verify that our new function works as expected. Same list L as before,")
println("There are $L_high \"high\" entries")
println("      and $L_low \"low\" entries")

println("Let's do the whole thing again. Again, running $n_trials trials on lists of length $n_rand")
unsorted_time = 0.0
sorted_time = 0.0
time_to_sort = 0.0
count_by_sort_time = 0.0
for trial = 1:n_trials
    R = rand(n_rand)
    tic()
    count_high_low(R)
    unsorted_time += toq()

    tic()
    count_high_low_by_sort(R)
    count_by_sort_time += toq()

    tic()
    sort!(R)
    time_to_sort += toq()

    # tic()
    # count_high_low(R)
    # sorted_time += toq()
end
println("Approximate total run time: $(unsorted_time + count_by_sort_time + time_to_sort) s")

println()
println("Average time to run count_high_low:               " * (@sprintf "%0.10f s" (unsorted_time/n_trials)) * " (a)")
println("Average time to run count_high_low_by_sort:       " * (@sprintf "%0.10f s" (count_by_sort_time/n_trials)) * " (b)")
println("Average time to sort the list:                    " * (@sprintf "%0.10f s" (time_to_sort/n_trials)) * " (c)")
println("Average time to perform binary search:            " * (@sprintf "%0.10f s" ((count_by_sort_time - time_to_sort)/n_trials)) * " (b) - (c)")
println()

println("""
Again, the speedup is overshadowed by the time to sort the list. But,
subtracting out the average time it takes to sort the list, we get a speedup of
 (a)/((b) - (c)) =
Speedup from binary search:  *** $(unsorted_time/(count_by_sort_time - time_to_sort)) ***

Compare to the speedup we got just from using the naive version of the function, with a sorted list:

Speedup from sorting:        *** $(sort_speedup) ***
""")

println("""
What I've been talking about is really the asymptotic behvaior of these
algorithms. Change the value of n_rand to 10^5 on line $(rand_line) and see what
happens to all these numbers.
""")
