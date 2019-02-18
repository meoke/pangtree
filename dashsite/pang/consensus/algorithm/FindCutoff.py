from abc import ABC
from typing import List
import numpy as np


class FindCutoff(ABC):
    @staticmethod
    def check_if_not_empty(compatibilities):
        if not list(compatibilities):
            raise ValueError("Empty compatibilities list. Cannot find cutoff.")

    @staticmethod
    def sort_and_get_value_following_max_distance(values):
        if len(values) == 1:
            return values[0]
        sorted_values = sorted(values)
        distances = np.array([sorted_values[i + 1] - sorted_values[i] for i in range(len(sorted_values) - 1)])
        max_distance_index = np.argmax(distances)
        return sorted_values[max_distance_index + 1]


class FindMaxCutoff(FindCutoff):
    def find_max_cutoff(self, compatibilities: List[float]) -> float:
        pass


class FindNodeCutoff(FindCutoff):
    def find_node_cutoff(self, compatibilities: List[float], so_far_cutoffs: List[float]) -> float:
        pass

    @staticmethod
    def use_max2(compatibilities):
        strategy = MAX2()
        return strategy.find_max_cutoff(compatibilities)


class MAX1(FindMaxCutoff):
    def __init__(self, cutoff_search_range: List[float]):
        if len(cutoff_search_range) != 2:
            raise ValueError("Cutoff search range must have length 2.")
        elif cutoff_search_range[1] < cutoff_search_range[0]:
            raise ValueError("For cutoff search range [x, y] x must be <= y.")
        self.cutoff_search_range = cutoff_search_range

    def find_max_cutoff(self, compatibilities: List[float]) -> float:
        FindCutoff.check_if_not_empty(compatibilities)
        min_search_pos = round((len(compatibilities) - 1) * self.cutoff_search_range[0])
        max_search_pos = round((len(compatibilities) - 1) * self.cutoff_search_range[1])
        sorted_comp = sorted(compatibilities)
        if min_search_pos == max_search_pos:
            return sorted_comp[min_search_pos]

        search_range = sorted(set(sorted_comp[min_search_pos: max_search_pos + 1]))

        return FindCutoff.sort_and_get_value_following_max_distance(search_range)


class MAX2(FindMaxCutoff):
    def find_max_cutoff(self, compatibilities: List[float]) -> float:
        FindCutoff.check_if_not_empty(compatibilities)
        return FindCutoff.sort_and_get_value_following_max_distance(compatibilities)


class NODE1(FindNodeCutoff):
    def __init__(self, multiplier):
        self.multiplier = multiplier


    def find_node_cutoff(self, compatibilities: List[float], so_far_cutoffs: List[float]) -> float:
        FindCutoff.check_if_not_empty(compatibilities)

        if len(compatibilities) == 1:
            return compatibilities[0]

        sorted_comp = sorted(compatibilities)
        mean_distance = (sorted_comp[-1] - sorted_comp[0]) / (len(sorted_comp) - 1)
        required_gap = mean_distance * self.multiplier
        cutoff = NODE1.get_value_following_first_gap_greater_than_required_gap(sorted_comp, required_gap)
        if cutoff is None:
            cutoff = NODE1.get_value_following_first_gap_greater_than_required_gap(sorted_comp, mean_distance)
        return cutoff

    @staticmethod
    def get_value_following_first_gap_greater_than_required_gap(sorted_comp, required_gap):
        if len(sorted_comp) == 1:
            return sorted_comp[0]
        distances = np.array([sorted_comp[i + 1] - sorted_comp[i] for i in range(len(sorted_comp) - 1)])
        if any(distances >= required_gap):
            a = np.where(distances >= required_gap)[0][0] + 1
            return sorted_comp[a]
        return None


class NODE2(FindNodeCutoff):
    def __init__(self, multiplier):
        self.multiplier = multiplier

    def find_node_cutoff(self, compatibilities: List[float], so_far_cutoffs: List[float]) -> float:
        if not so_far_cutoffs:
            return self.use_node1(compatibilities)

        guard = min(so_far_cutoffs)
        sorted_comp = sorted(compatibilities)
        if guard <= sorted_comp[0]:
            return sorted_comp[0]
        elif guard > sorted_comp[-1]:
            return self.use_node1(compatibilities)
        else:
            sorted_comp = sorted(sorted_comp + [guard])
            mean_distance = (sorted_comp[-1] - sorted_comp[0]) / (len(sorted_comp) - 1)
            required_gap = mean_distance * self.multiplier
            if sorted_comp.count(guard) > 1:
                cutoff_search_area = [c for c in sorted_comp if c <= guard]
            else:
                cutoff_search_area = [c for c in sorted_comp if c < guard]
            cutoff = NODE1.get_value_following_first_gap_greater_than_required_gap(cutoff_search_area, required_gap)
            if cutoff is None:
                comp_greater_than_guard = [c for c in sorted_comp if c > guard]
                cutoff = min(comp_greater_than_guard)
            return cutoff

    def use_node1(self, compatibilities):
        strategy = NODE1(self.multiplier)
        return strategy.find_node_cutoff(compatibilities, [])


class NODE3(FindNodeCutoff):
    def find_node_cutoff(self, compatibilities: List[float], so_far_cutoffs: List[float]) -> float:
        if not so_far_cutoffs:
            return self.use_max2(compatibilities)

        guard = min(so_far_cutoffs)
        sorted_comp = sorted(compatibilities)
        if guard <= sorted_comp[0]:
            return sorted_comp[0]
        elif guard >= sorted_comp[-1]:
            return self.use_max2(compatibilities)
        else:
            first_comp_greater_than_guard_index = [i for i, c in enumerate(sorted_comp) if c > guard][0]
            return self.use_max2(sorted_comp[0:first_comp_greater_than_guard_index + 1])


class NODE4(FindNodeCutoff):
    def find_node_cutoff(self, compatibilities: List[float], so_far_cutoffs: List[float]):
        return self.use_max2(compatibilities)


