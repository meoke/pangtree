from .Direction import Direction


class Block:
    def __init__(self, block_id: int, order_id: int, direction: Direction, alignment):
        #TODO uzupełnić typ dla argumentu alignment
        self.block_id = block_id
        self.order_id = order_id
        self.direction = direction
        self.alignment = alignment
