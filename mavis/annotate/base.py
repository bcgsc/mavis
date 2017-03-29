from ..interval import Interval
import re
import itertools


class ReferenceName(str):
    def __eq__(self, other):
        putative_other = [other, 'chr' + str(other), re.sub('^chr', '', str(other))]
        putative_self = [self, 'chr' + str(self), re.sub('^chr', '', str(self))]

        for s, o in itertools.product(putative_self, putative_other):
            if str.__eq__(s, o):
                return True
        return False

    def __hash__(self):
        return hash(re.sub('^chr', '', str(self)))


class BioInterval:
    def __init__(self, reference_object, start, end=None, name=None, seq=None, data=None):
        """
        Args:
            reference_object: the object this interval is on
            start (int) start of the interval (inclusive)
            end (int): end of the interval (inclusive)
            name: optional
            seq (str): the seq relating to this interval

        Example:
            >>> b = BioInterval('1', 12572784, 12578898, 'q22.2')
            >>> b[0]
            12572784
            >>> b[1]
            12578898
        """
        start = int(start)
        end = int(end) if end is not None else None
        data = {} if data is None else data
        self.reference_object = reference_object
        self.name = name
        self.position = Interval(start, end, number_type=int)
        self.seq = seq if not seq else seq.upper()
        self.data = {}
        self.data.update(data)

    @property
    def start(self):
        """*int*: the start position"""
        return self.position.start

    @property
    def end(self):
        """*int*: the end position"""
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
        return self.position.length()

    def key(self):
        """:class:`tuple`: a tuple representing the items expected to be unique. for hashing and comparing"""
        return (self.reference_object, self.position, self.seq, self.name)

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

    def get_seq(self, REFERENCE_GENOME=None, ignore_cache=False):
        raise NotImplementedError('abstract method must be overidden')

    def get_strand(self):
        """
        pulls strand information from the current object, or follows reference
        objects until the strand is found

        Returns:
            STRAND: the strand of this or any of its reference objects

        Raises:
            AttributeError: raised if the strand is not set on this or any of its reference objects
        """
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
        raise AttributeError('strand has not been defined', self)

    def get_chr(self):
        """
        pulls chromosome information from the current object, or follows reference
        objects until the chromosome is found

        Returns:
            str: the chromosome of this or any of its reference objects

        Raises:
            AttributeError: raised if the chromosome is not set on this or any of its reference objects
        """
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

    def to_dict(self):
        """
        creates a dictionary representing the current object

        Returns:
            :class:`dict` by :class:`str`: the dictionary of attribute values
        """
        d = {
            'name': self.name,
            'start': self.start,
            'end': self.end,
            'type': self.__class__.__name__,
            'seq': self.seq,
            'data': self.data
        }
        try:
            d['reference_object'] = self.reference_object.name
        except AttributeError:
            d['reference_object'] = str(self.reference_object)
        return d

    def __repr__(self):
        cls = self.__class__.__name__
        refname = self.reference_object
        try:
            refname = self.reference_object.name
        except AttributeError:
            pass
        return '{}({}:{}-{}, name={})'.format(cls, refname, self.start, self.end, self.name)
