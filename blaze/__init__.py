from .blaze import main as _pipeline

def blaze(argv = None):
    if argv == None:
        _pipeline('-h')
        return 
    else: 
        return _pipeline(argv)

