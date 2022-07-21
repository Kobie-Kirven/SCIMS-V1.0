# Holds error classes

class Error(Exception):
    """Base class for other exceptions"""
    pass

class WrongFile(Error):
    """Raised when the SAM or BAM file doesn't look right"""
    pass

class DoesNotExist(Error):
    """
    Raised when the input file does not exist
    """
    pass