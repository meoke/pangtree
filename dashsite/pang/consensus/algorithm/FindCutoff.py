import csv
from abc import ABC
from typing import List
import numpy as np


class FindCutoff(ABC):
    def __init__(self, cutoffs_log_file_path=None):
        self.cutoffs_log_file_path = cutoffs_log_file_path
        self.log_counter = -1
        if self.cutoffs_log_file_path:
            with open(self.cutoffs_log_file_path, 'a') as output:
                csv_writer = csv.writer(output, delimiter=',')
                csv_writer.writerow(["id",
                                     "strategy",
                                     "compatibilities",
                                     "params",
                                     "cutoff",
                                     "reason"
                                     ])

    def log_cutoff_search(self, compatibilities, cutoff, reason, params):
        self.log_counter += 1
        with open(self.cutoffs_log_file_path, 'a') as output:
            csv_writer = csv.writer(output, delimiter=',')
            csv_writer.writerow([self.log_counter,
                                 self.__class__.__name__,
                                 compatibilities,
                                 params,
                                 cutoff,
                                 reason
                                 ])

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
    def find_max_cutoff(self, compatibilities: List[float], log: bool=False) -> float:
        pass


class FindNodeCutoff(FindCutoff):
    def find_node_cutoff(self, compatibilities: List[float], so_far_cutoffs: List[float], log: bool = False) -> float:
        pass

    @staticmethod
    def get_max2_result(compatibilities):
        strategy = MAX2()
        return strategy.find_max_cutoff(compatibilities, log=False)


class MAX1(FindMaxCutoff):
    def __init__(self, cutoff_search_range: List[float], cutoffs_log_file_path = None):
        super().__init__(cutoffs_log_file_path)
        if len(cutoff_search_range) != 2:
            raise ValueError("Cutoff search range must have length 2.")
        elif cutoff_search_range[1] < cutoff_search_range[0]:
            raise ValueError("For cutoff search range [x, y] x must be <= y.")
        self.cutoff_search_range = cutoff_search_range

    def find_max_cutoff(self, compatibilities: List[float], log: bool=False) -> float:
        FindCutoff.check_if_not_empty(compatibilities)
        min_search_pos = round((len(compatibilities) - 1) * self.cutoff_search_range[0])
        max_search_pos = round((len(compatibilities) - 1) * self.cutoff_search_range[1])
        sorted_comp = sorted(compatibilities)
        if min_search_pos == max_search_pos:
            cutoff = sorted_comp[min_search_pos]
            reason = "Single value remained after applying cutoff search range."
        else:
            search_range = sorted(sorted_comp[min_search_pos: max_search_pos + 1])
            cutoff = FindCutoff.sort_and_get_value_following_max_distance(search_range)
            reason = "Value after max distance in cutoff search range."

        if log:
            self.log_cutoff_search(compatibilities, cutoff, reason, f"r: {self.cutoff_search_range}")
        return cutoff


class MAX2(FindMaxCutoff):
    def find_max_cutoff(self, compatibilities: List[float], log: bool=False) -> float:
        FindCutoff.check_if_not_empty(compatibilities)
        cutoff = FindCutoff.sort_and_get_value_following_max_distance(compatibilities)
        reason = "Value after max distance in cutoff search range."
        if log:
            self.log_cutoff_search(compatibilities, cutoff, reason, "")
        return cutoff


class NODE1(FindNodeCutoff):
    def __init__(self, multiplier, cutoffs_log_file_path = None):
        super().__init__(cutoffs_log_file_path)
        self.multiplier = multiplier

    def find_node_cutoff(self, compatibilities: List[float], so_far_cutoffs: List[float], log: bool=False) -> float:
        FindCutoff.check_if_not_empty(compatibilities)

        if len(compatibilities) == 1:
            cutoff = compatibilities[0]
            reason = "Single compatibility"
        else:
            sorted_comp = sorted(compatibilities)
            # sorted_comp = sorted(set(compatibilities))
            mean_distance = (sorted_comp[-1] - sorted_comp[0]) / (len(sorted_comp) - 1)
            required_gap = mean_distance * self.multiplier
            cutoff = NODE1.get_value_following_first_gap_greater_than_required_gap(sorted_comp, required_gap)
            reason = "Value after mean_distance * multiplier"
            if cutoff is None:
                cutoff = NODE1.get_value_following_first_gap_greater_than_required_gap(sorted_comp, mean_distance)
                reason = "Value after mean_distance * 1"
        if log:
            self.log_cutoff_search(compatibilities, cutoff, reason, f"multiplier: {self.multiplier}")
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
    def __init__(self, multiplier, cutoffs_log_file_path = None):
        super().__init__(cutoffs_log_file_path)
        self.multiplier = multiplier

    def find_node_cutoff(self, compatibilities: List[float], so_far_cutoffs: List[float], log: bool=False) -> float:
        if not so_far_cutoffs:
            cutoff = self.get_node1_result(compatibilities)
            reason = "No so far cutoffs. First child in node."
        else:
            guard = min(so_far_cutoffs)
            sorted_comp = sorted(compatibilities)
            if guard <= sorted_comp[0]:
                cutoff = sorted_comp[0]
                reason = "guard <= min(compatibilities), return min(compatibilities)"
            elif guard > sorted_comp[-1]:
                cutoff = self.get_node1_result(compatibilities)
                reason = "guard > max(compatibilities), use node 1."
            else:
                sorted_comp = sorted(sorted_comp + [guard])
                mean_distance = (sorted_comp[-1] - sorted_comp[0]) / (len(sorted_comp) - 1)
                required_gap = mean_distance * self.multiplier
                if sorted_comp.count(guard) > 1:
                    cutoff_search_area = [c for c in sorted_comp if c <= guard]
                else:
                    cutoff_search_area = [c for c in sorted_comp if c < guard]
                cutoff = NODE1.get_value_following_first_gap_greater_than_required_gap(cutoff_search_area, required_gap)
                reason = "Value after mean_distance * multiplier"
                if cutoff is None:
                    comp_greater_than_guard = [c for c in sorted_comp if c > guard]
                    cutoff = min(comp_greater_than_guard)
                    reason = "No gap greater then mean_distance*multiplier. Take first to the right."
        if log:
            self.log_cutoff_search(compatibilities, cutoff, reason, f"multiplier: {self.multiplier}, so_far_cutoffs: {so_far_cutoffs}")
        return cutoff

    def get_node1_result(self, compatibilities):
        strategy = NODE1(self.multiplier)
        return strategy.find_node_cutoff(compatibilities, [], log=False)


class NODE3(FindNodeCutoff):
    def find_node_cutoff(self, compatibilities: List[float], so_far_cutoffs: List[float], log: bool=False) -> float:
        if not so_far_cutoffs:
            cutoff = self.get_max2_result(compatibilities)
            reason = "No so far cutoffs. Use max 2."
        else:
            guard = min(so_far_cutoffs)
            sorted_comp = sorted(compatibilities)
            if guard <= sorted_comp[0]:
                cutoff = sorted_comp[0]
                reason = "guard < min(compatibilities). Return min(compatibilities)."
            elif guard >= sorted_comp[-1]:
                cutoff = self.get_max2_result(compatibilities)
                reason = "guard > max(compatibilities). Use max 2."
            else:
                first_comp_greater_than_guard_index = [i for i, c in enumerate(sorted_comp) if c > guard][0]
                cutoff = self.get_max2_result(sorted_comp[0:first_comp_greater_than_guard_index + 1])
                reason = "Use max 2 on sorted_comp[0:first_comp_greater_than_guard_index + 1]"
        if log:
            self.log_cutoff_search(compatibilities, cutoff, reason, f"so far cutoffs: {so_far_cutoffs}")
        return cutoff


class NODE4(FindNodeCutoff):
    def find_node_cutoff(self, compatibilities: List[float], so_far_cutoffs: List[float], log=False):
        cutoff = self.get_max2_result(compatibilities)
        reason = "Use MAX2"
        if log:
            self.log_cutoff_search(compatibilities, cutoff, reason, "")
        return cutoff
