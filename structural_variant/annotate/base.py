from ..interval import Interval


class BioInterval:
    def __init__(self, reference_object, start, end, name=None):
        """
        Args:
            reference_object: the object this interval is on
            start (int) start of the interval (inclusive)
            end (int): end of the interval (inclusive)
            name: optional

        Example:
            >>> b = BioInterval('1', 12572784, 12578898, 'q22.2')
            >>> b[0]
            12572784
            >>> b[1]
            12578898
        """
        self.reference_object = reference_object
        self.name = name
        self.position = Interval(start, end)

    @property
    def start(self):
        return self.position.start

    @property
    def end(self):
        return self.position.end

    def __getitem__(self, index):
        return Interval.__getitem__(self, index)

    def __len__(self):
        """
        Example:
            >>> b = BioInterval('1', 12572784, 12578898, 'q22.2')
            >>> len(b)
            6115
        """
        return len(self.position)

    @property
    def key(self):
        return self.reference_object, self.position, self.name

    def __eq__(self, other):
        if not hasattr(other, 'key'):
            return False
        else:
            return self.key == other.key

    def __lt__(self, other):
        if other.reference_object != self.reference_object:
            if self.reference_object < other.reference_object:
                return True
            else:
                return False
        elif self.position < other.position:
            return True
        else:
            return False

    def __hash__(self):
        return hash(self.key)

    def get_sequence(self, REFERENCE_GENOME=None):
        raise NotImplementedError('abstract method must be overidden')
