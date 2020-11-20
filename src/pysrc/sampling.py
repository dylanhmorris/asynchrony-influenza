import random
from collections import Counter


def weighted_choice(weights):
    """
    Fast random choice using non-normalized weights
    implemented using a python algorithm
    developed by Eli Bendersky
    (http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/)
    """
    rnd = random.random() * sum(weights)
    for index, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return index


def _cached_weighted_choice(weights, weight_sum):
    """
    Random choice using non-normalized weights
    implemented using a python algorithm
    developed by Eli Bendersky
    (http://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/)

    Optimized for repeated use with the same weight set in that it
    it takes a precalculated weight_sum alongside the weights
    instead of calculating it from the weights at runtime
    """
    rnd = random.random() * weight_sum
    for index, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return index


def weighted_draw(sample_space_counts, size):
    sorted_counts = sample_space_counts.most_common()
    weights = [item[1] for item in sorted_counts]
    output = Counter({item[0]: 0 for item in sorted_counts})
    weight_sum = sum(weights)
    if weight_sum == 0:
        print("WARNING: No population to sample from. "
              "Returning an empty Counter. ")
        return output
    for k in range(size):
        choice = _cached_weighted_choice(weights, weight_sum)
        strain = sorted_counts[choice][0]
        weights[choice] += -1
        weight_sum += - 1
        output[strain] += 1
        if weight_sum == 0:
            print("WARNING: Sampled population down to zero.")
            break

    return output


def weighted_draw_pure_python(sample_space_counts, size):
    sorted_counts = sample_space_counts.most_common()
    weights = [item[1] for item in sorted_counts]
    output = Counter({item[0]: 0 for item in sorted_counts})
    weight_sum = sum(weights)
    if weight_sum == 0:
        print("WARNING: No population to sample from. "
              "Returning an empty Counter. ")
        return output
    for k in range(size):
        choice = _cached_weighted_choice(weights, weight_sum)
        strain = sorted_counts[choice][0]
        weights[choice] += -1
        weight_sum += - 1
        output[strain] += 1
        if weight_sum == 0:
            print("WARNING: Sampled population down to zero.")
            break

    return output


def big_set_sample(set_to_sample, set_size=None):
    """
    This was a nice idea, but it fails because
    set iteration ordering is neither perfectly random
    nor perfectly deterministic for a given state of the
    set
    """

    if not set_size:
        set_size = len(set_to_sample)

    generator = (item for item in set_to_sample)
    mass = 1
    increment = 1/set_size
    for iteration, value in enumerate(generator):
        mass -= increment
        if random.random() > mass:
            return value
    raise BaseException("Failed to return a value. Check set size")
