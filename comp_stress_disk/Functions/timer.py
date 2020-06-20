import time
from functools import wraps

def timer(func):
    """A decorator that prints how long a function took to run."""
    # Define the wrapper function to return.
    @wraps(func)
    def wrapper(*args, **kwargs):
        # When wrapper() is called, get the current time.
        t_start = time.time()
        # Call the decorated function and store the result.
        result = func(*args, **kwargs)
        # Get the total time it took to run, and print it.
        t_total = time.time() - t_start
        print('{} took {}s'.format(func.__name__, t_total))
    return result
return wrapper