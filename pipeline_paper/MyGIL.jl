
module MyGIL

const PYLOCK = Ref{ReentrantLock}()

function __init__() # instantiate the lock
    PYLOCK[] = ReentrantLock()
end

# acquire the lock before any code calls Python
pylock(f::Function) = Base.lock(PYLOCK[]) do
    prev_gc = GC.enable(false)
    try
        return f()
    finally
        GC.enable(prev_gc) # recover previous state
    end
end

end # m