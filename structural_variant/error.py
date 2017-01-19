

class DiscontinuousMappingError(Exception):
    """
    class to describe errors associated with translating
    coordinates between continuous and discontinuous intervals
    """

    def __init__(self, *pos, before=None, after=None):
        self.before = before
        self.after = after
        if self.before is None and self.after is None:
            raise AttributeError('please specify before OR after OR both')

        self.msg = ' '.join([str(x) for x in pos])

    def __str__(self):
        name = self.__class__.__name__
        return '{0}<[{1},{2}], {3}>'.format(name, self.before, self.after, self.msg)


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

    def __init__(self, *pos):
        self.msg = ' '.join([str(p) for p in pos])

    def __str__(self):
        name = self.__class__.__name__
        return '{0}<rearrangement would not produce a proper genetic molecule: {1}>'.format(name, self.msg)
