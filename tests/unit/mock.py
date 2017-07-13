import types


class Mock:
    def __init__(self, **kwargs):
        for attr, val in kwargs.items():
            setattr(self, attr, val)

    def bind_method(self, **kwargs):
        for attr, val in kwargs.items():
            val = types.MethodType(val, self)  # bind the method to self
            setattr(self, attr, val)

    def add_attr(self, attr, val):
        setattr(self, attr, val)

    def __contains__(self, item):
        if hasattr(self, item):
            return True
        return False


class MockFunction:
    def __init__(self, return_value):
        self.return_value = return_value

    def __call__(self, *pos, **kwargs):
        return self.return_value


class MockLongString:
    def __init__(self, string, offset):
        self.string = string
        self.offset = offset

    def __len__(self):
        return len(self.string) + self.offset

    def __getitem__(self, index):
        if not isinstance(index, slice):
            index = slice(index, index + 1)
        index = slice(
            index.start - self.offset,
            index.stop - self.offset,
            index.step)
        if index.start < 0:
            raise NotImplementedError('string portion not given')
        return self.string[index]
