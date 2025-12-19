#
# A sort of progress bar with time stamps
#

struct Progress
    niterations::Integer
    divisor::Float64
    interval::Integer
    start_time
    function Progress(niterations; divisor=10)
        interval = ceil(Int, niterations/divisor)
        new(niterations, divisor, interval, now())
    end
end

function print_progress(progress::Progress, count)
    printed = false
    if(count %  progress.interval == 0 || count == progress.niterations)
        elapsed_time = Dates.value(now()-progress.start_time)
        hours = floor(Int, elapsed_time/(1000*3600))
        minutes = floor(Int, elapsed_time/(1000*60)) % 60
        seconds = floor(Int, elapsed_time/1000) % 60
        timestr = string(hours) * ":" * string(minutes) * ":" * string(seconds) 
        println(round(100*count/progress.niterations), " % at ",timestr)
        flush(stdout)
        printed = true
    end
    return printed
end
