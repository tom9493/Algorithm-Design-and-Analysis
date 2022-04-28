import random

def prime_test(N, k):
    # This is main function, that is connected to the Test button. You don't need to touch it.
    return fermat(N, k), miller_rabin(N, k)

#a , N-1 , N
def mod_exp(x, y, N):                       # Whole function O(n^3)
    if y == 0:
        return 1
    z = mod_exp(x, y//2, N)                 # O(y) or O(n-1) or O(n) for how many calls, O(n^2) for division
    if y % 2 == 0:                          # O(n^2) divide
        return z**2%N                       # O(n^2) multiply
    else:
        return x*z**2%N                     # O(n^2) multiply


def fprobability(k):
    return 1 - 2**-k                        # Minus is Big O(n), multiply is O(n^2)


def mprobability(k):
    return 1 - 4**-k                        # O(n^2) multiply


def fermat(N, k):                           # Whole thing O(n^3)
    for i in range(k):                      # O(k)
        a = random.randint(1, N - 1)        # O(n)
        f = mod_exp(a, N - 1, N)            # O(n^3)
        if f != 1:
            return 'composite'
    return 'prime'


def miller_rabin(N, k):                      # Whole function is O(n^3)
    e = N - 1                                # O(n) addition
    if e % 2 != 0:                           # O(n^2) divide
        return 'composite'
    for i in range(k):                       # O(k) for loop
        a = random.randint(2, N - 2)         # O(n) subtract
        while e % 2 == 0:                    # O(n^2) divide
            f = mod_exp(a, e, N)             # O(n^2)
            if f == -1 or f == N - 1:        # O(2n) 2 subtractions
                continue
            if f == 1:
                e //= 2                      # O(n^2) divide
            else:
                return 'composite'
    return 'prime'
