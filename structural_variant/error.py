
class DiscontiuousMappingError(Exception):
    """ 
    class to describe errors associated with translating 
    coordinates between continuous and discontinuous intervals
    """
    def __init__(self, msg, *pos, **kwargs):
        self.between = kwargs.pop('between', False)
        self.before = kwargs.pop('before', False)
        self.after = kwargs.pop('after', False)

        if self.between is not False:
            self.pos = self.between
            self.between = True
        elif self.after is not False:
            self.pos = self.after
            self.after = True
        elif self.before is not False:
            self.pos = self.before
            self.before = True

        if sum([ 1 for arg in [self.between, self.after, self.before] if arg]) != 1:
            raise AttributeError('before, after, and between arguments '
                    'are both required and mutually exclusive')
        self.msg = ' '.join([msg] + [p for p in pos])
        if kwargs:
            raise AttributeError('unexpected keyword argument', kwargs)

    def __str__(self):
        name = self.__class__.__name__
        if self.between:
            return '{0}<between={1}, {2}>'.format(name, self.pos, self.msg)
        elif self.before:
            return '{0}<before={1}, {2}>'.format(name, self.pos, self.msg)
        elif self.after:
            return '{0}<after={1}, {2}>'.format(name, self.pos, self.msg)

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
        self.msg = ' '.join(list(pos))

    def __str__(self):
        name = self.__class__.__name__
        return '{0}<rearrangement would not produce a proper genetic molecule: {1}>'.format(name, self.msg) 
