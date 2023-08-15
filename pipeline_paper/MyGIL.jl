
module MyGIL

const PYLOCK = Ref{ReentrantLock}()

function __init__() # instantiate the lock
    PYLOCK[] = ReentrantLock()
end

# acquire the lock before any code calls Python
pylock(f::Function) = Base.lock(f, PYLOCK[])


end # m