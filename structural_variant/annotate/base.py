from ..interval import Interval


class BioInterval:
    def __init__(self, reference_object, start, end, name=None, sequence=None):
        """
        Args:
            reference_object: the object this interval is on
            start (int) start of the interval (inclusive)
            end (int): end of the interval (inclusive)
            name: optional
            sequence (str): the sequence relating to this interval

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
        self.sequence = sequence if not sequence else sequence.upper()

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

    def key(self):
        return (self.reference_object, self.position, self.sequence, self.name)

    def __eq__(self, other):
        if not hasattr(other, 'key'):
            return False
        else:
            return self.key() == other.key()

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
        return hash(self.key())

    def get_sequence(self, REFERENCE_GENOME=None, ignore_cache=False):
        raise NotImplementedError('abstract method must be overidden')

    def get_strand(self):
        try:
            if self.strand is not None:
                return self.strand
        except AttributeError:
            pass
        tried = set()
        parent = self.reference_object
        while True:
            if parent in tried:
                break
            try:
                if parent.strand is not None:
                    return parent.strand
            except AttributeError:
                pass
            tried.add(parent)
            try:
                parent = parent.reference_object
            except AttributeError:
                break
        raise AttributeError('strand has not been defined')
    
    def get_chr(self):
        try:
            if self.chr is not None:
                return self.chr
        except AttributeError:
            pass

        tried = set()
        parent = self.reference_object
        while True:
            if parent in tried:
                break
            try:
                if parent.chr is not None:
                    return parent.chr
            except AttributeError:
                pass
            tried.add(parent)
            try:
                parent = parent.reference_object
            except AttributeError:
                break
        raise AttributeError('chr has not been defined')
