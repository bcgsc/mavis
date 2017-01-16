

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


class StrandSpecificityError(Exception):
    """
    raised when STRAND.NS is used a give process requires that the strand be specified
    """

    def __init__(self, *pos):
        self.msg = ' '.join(list(pos))

    def __str__(self):
        name = self.__class__.__name__
        return '{0}<strand must be specified: {1}>'.format(name, self.msg)


class InvalidRearrangement(Exception):

    def __init__(self, *pos):
        self.msg = ' '.join([str(p) for p in pos])

    def __str__(self):
        name = self.__class__.__name__
        return '{0}<rearrangement would not produce a proper genetic molecule: {1}>'.format(name, self.msg)
