
module MyGIL

const PYLOCK = Ref{ReentrantLock}()

function __init__() # instantiate the lock
    PYLOCK[] = ReentrantLock()
end

# acquire the lock before any code calls Python
# pylock(f::Function) = Base.lock(f, PYLOCK[])
function pylock(f::Function)
    GC.enable(false)
    lock(PYLOCK[])
    try
        return f()
    finally
        unlock(PYLOCK[])
        GC.enable(true)
    end
end

end # m