

class NotSpecifiedError(Exception):
    """
    raised when information is required for a function but has not been given

    for example if strand was required but had been set to STRAND.NS then this
    error would be raised
    """
    pass


class DrawingFitError(Exception):
    pass


class InvalidRearrangement(Exception):
    pass
