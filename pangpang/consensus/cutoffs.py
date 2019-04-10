from abc import ABC
from typing import List, Optional
import numpy as np

from pangpang.consensus.ConsensusTree import CompatibilityToPath
from pangpang.consensus import input_types


class FindCutoffResult:
    def __init__(self, cutoff: CompatibilityToPath, explanation: str):
        self.cutoff: CompatibilityToPath = cutoff
        self.explanation: str = explanation


class FindCutoff(ABC):
    @staticmethod
    def break_if_empty(compatibilities: List[CompatibilityToPath]) -> None:
        if not list(compatibilities):
            raise ValueError("Empty compatibilities list. Cannot find cutoff.")

    @staticmethod
    def sort_and_get_value_following_max_distance(values: List[CompatibilityToPath]) -> CompatibilityToPath:
        if len(values) == 1:
            return values[0]
        sorted_values = sorted(values)
        distances = np.array([sorted_values[i + 1] - sorted_values[i] for i in range(len(sorted_values) - 1)])
        max_distance_index: int = np.argmax(distances)
        return sorted_values[max_distance_index + 1]


class FindMaxCutoff(FindCutoff):
    def find_max_cutoff(self, compatibilities: List[CompatibilityToPath]) -> FindCutoffResult:
        pass


class FindNodeCutoff(FindCutoff):
    def find_node_cutoff(self,
                         compatibilities: List[CompatibilityToPath],
                         so_far_cutoffs: List[CompatibilityToPath]) -> FindCutoffResult:
        pass

    @staticmethod
    def get_max2_result(compatibilities: List[CompatibilityToPath]) -> CompatibilityToPath:
        strategy = MAX2()
        max2_result = strategy.find_max_cutoff(compatibilities)
        return max2_result.cutoff


class MAX1(FindMaxCutoff):
    def __init__(self, cutoff_search_range: input_types.Range):
        super().__init__()
        # todo remove below data check
        if len(cutoff_search_range.value) != 2:
            raise ValueError("Cutoff search range must have length 2.")
        elif cutoff_search_range.value[1] < cutoff_search_range.value[0]:
            raise ValueError("For cutoff search range [x, y] x must be <= y.")
        self.cutoff_search_range: input_types.Range = cutoff_search_range

    def find_max_cutoff(self, compatibilities: List[CompatibilityToPath]) -> FindCutoffResult:
        FindCutoff.break_if_empty(compatibilities)
        min_search_pos = round((len(compatibilities) - 1) * self.cutoff_search_range.value[0])
        max_search_pos = round((len(compatibilities) - 1) * self.cutoff_search_range.value[1])
        sorted_comp = sorted(compatibilities)
        if min_search_pos == max_search_pos:
            cutoff = sorted_comp[min_search_pos]
            reason = "Single value remained after applying cutoff search range."
        else:
            search_range = sorted(sorted_comp[min_search_pos: max_search_pos + 1])
            cutoff = FindCutoff.sort_and_get_value_following_max_distance(search_range)
            reason = "Value after max distance in cutoff search range."

        return FindCutoffResult(cutoff, reason)


class MAX2(FindMaxCutoff):
    def find_max_cutoff(self, compatibilities: List[CompatibilityToPath]) -> FindCutoffResult:
        FindCutoff.break_if_empty(compatibilities)
        cutoff = FindCutoff.sort_and_get_value_following_max_distance(compatibilities)
        reason = "Value after max distance in cutoff search range."
        return FindCutoffResult(cutoff, reason)


class NODE1(FindNodeCutoff):
    def __init__(self, multiplier: input_types.Multiplier):
        super().__init__()
        self.multiplier: input_types.Multiplier = multiplier

    def find_node_cutoff(self, compatibilities: List[CompatibilityToPath],
                         so_far_cutoffs: List[CompatibilityToPath]) -> FindCutoffResult:
        FindCutoff.break_if_empty(compatibilities)

        if len(compatibilities) == 1:
            cutoff = compatibilities[0]
            reason = "Single compatibility"
        else:
            sorted_comp = sorted(compatibilities)
            mean_distance = (sorted_comp[-1] - sorted_comp[0]).value / (len(sorted_comp) - 1)
            required_gap = mean_distance * self.multiplier.value
            cutoff = NODE1.get_value_following_first_gap_greater_than_required_gap(sorted_comp, required_gap)
            reason = "Value after mean_distance * multiplier"
            if cutoff is None:
                cutoff = NODE1.get_value_following_first_gap_greater_than_required_gap(sorted_comp, mean_distance)
                reason = "Value after mean_distance * 1"
        return FindCutoffResult(cutoff, reason)

    @staticmethod
    def get_value_following_first_gap_greater_than_required_gap(sorted_comp: List[CompatibilityToPath],
                                                                required_gap: float) -> Optional[CompatibilityToPath]:
        if len(sorted_comp) == 1:
            return sorted_comp[0]
        distances = np.array([(sorted_comp[i + 1] - sorted_comp[i]).value for i in range(len(sorted_comp) - 1)])
        if any(distances >= required_gap):
            a = np.where(distances >= required_gap)[0][0] + 1
            return sorted_comp[a]
        return None


class NODE2(FindNodeCutoff):
    def __init__(self, multiplier: input_types.Multiplier):
        super().__init__()
        self.multiplier = multiplier

    def find_node_cutoff(self,
                         compatibilities: List[CompatibilityToPath],
                         so_far_cutoffs: List[CompatibilityToPath]) -> FindCutoffResult:
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
                mean_distance = (sorted_comp[-1] - sorted_comp[0]).value / (len(sorted_comp) - 1)
                required_gap = mean_distance * self.multiplier.value
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
        return FindCutoffResult(cutoff, reason)

    def get_node1_result(self, compatibilities: List[CompatibilityToPath]) -> CompatibilityToPath:
        strategy = NODE1(self.multiplier)
        node1_result = strategy.find_node_cutoff(compatibilities, [])
        return node1_result.cutoff


class NODE3(FindNodeCutoff):
    def find_node_cutoff(self,
                         compatibilities: List[CompatibilityToPath],
                         so_far_cutoffs: List[CompatibilityToPath]) -> FindCutoffResult:
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
        return FindCutoffResult(cutoff, reason)


class NODE4(FindNodeCutoff):
    def find_node_cutoff(self,
                         compatibilities: List[CompatibilityToPath],
                         so_far_cutoffs: List[CompatibilityToPath]) -> FindCutoffResult:
        cutoff = self.get_max2_result(compatibilities)
        reason = "Use MAX2"
        return FindCutoffResult(cutoff, reason)
