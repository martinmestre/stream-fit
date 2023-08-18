
module MyGIL

const PYLOCK = Ref{ReentrantLock}()

function __init__() # instantiate the lock
    PYLOCK[] = ReentrantLock()
end

# acquire the lock before any code calls Python
<<<<<<< HEAD
# pylock(f::Function) = Base.lock(f, PYLOCK[])
function pylock(f::Function)
    GC.enable(false)
    lock(PYLOCK[])
    try
        return f()
    finally
        unlock(PYLOCK[])
        GC.enable(true)
=======
pylock(f::Function) = Base.lock(PYLOCK[]) do
    prev_gc = GC.enable(false)
    try
        return f()
    finally
        GC.enable(prev_gc) # recover previous state
>>>>>>> 1783fa9ce18502164e93bee370b7d5df8be48a38
    end
end

end # m