from typing import List


class FilePart(object):
    def __init__(self, content: List[str], lines_count: int):
        self.content = content
        self.lines_count = lines_count
