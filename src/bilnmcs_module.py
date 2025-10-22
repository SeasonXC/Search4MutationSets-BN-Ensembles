# lnmcs_module.py
import random
import time

# ------------------ Helper Functions ------------------ #

def normalize_sorted_list(state_list):
    """Return a NEW sorted list of (gene, bool) pairs."""
    return sorted(state_list, key=lambda x: (x[0], x[1]))

def normalize_key(state_list):
    """Hashable, order-independent cache key."""
    return tuple(normalize_sorted_list(state_list))

def is_terminal(state_set, depth):
    """Check if the mutation set has reached the desired size."""
    return len(state_set) >= depth

def legal_moves_fn(state_set, all_moves):
    """Return list of mutations not yet applied. Preserves order of all_moves."""
    applied_genes = {gene for gene, _ in state_set}
    return [m for m in all_moves if m[0] not in applied_genes]

def random_playout(state_set, all_moves, depth, ec):
    """Complete the mutation set randomly up to `depth` and score it."""
    remaining = legal_moves_fn(state_set, all_moves)
    k = depth - len(state_set)
    # Guard (should be safe if depth ≤ total unique genes)
    k = max(0, min(k, len(remaining)))
    tail = random.sample(remaining, k=k)
    full_set = normalize_sorted_list(list(state_set) + tail)
    score = ec.evaluate(full_set)
    return score, full_set
# ------------------ LNMCS (Lazy NMCS) ------------------ #
'''
def bilnmcs(state, level, depth, all_moves, ec_main,ec_fast, *,
          b=2, r=0.5, e=None, timeout_sec=None):
    """
    Lazy NMCS entrypoint.

    Args:
        state: current mutation set (list[(gene, bool)])
        level: nesting level (>= 0)
        depth: target set size (terminal when len(state) >= depth)
        all_moves: full ordered move list (same as best_moves in NMCS)
        ec: evaluator with ec.evaluate(full_set)
        b: number of cheap playouts per candidate move (small: 2–5)
        r: pruning ratio in [0,1]; higher = less pruning
        e: optional cap on number of moves to evaluate per node
        timeout_sec: optional wallclock limit

    Returns:
        (best_score, best_state_sorted)
    """
    deadline = time.time() + timeout_sec if timeout_sec else None

    # Per-depth running stats used for thresholds
    # tr[d] = {'mean': float, 'count': int}
    # trmax[d] = float (best mean seen at depth d)
    max_depth = depth  # depth index won’t exceed this
    tr     = [{'mean': 0.0, 'count': 0} for _ in range(max_depth + 1)]
    trmax  = [float('-inf') for _ in range(max_depth + 1)]

    best = {'score': -float('inf'), 'state': []}
    score, s = _bilnmcs(state, level, depth, all_moves, ec_main, ec_fast,
                      cache={}, deadline=deadline, best=best,
                      tr=tr, trmax=trmax, b=b, r=r, e=e)
    return score, normalize_sorted_list(s)


def _bilnmcs(state, level, depth, all_moves, ec_main, ec_fast, *,
           cache, deadline, best, tr, trmax, b, r, e):
    # Timeout guard
    if deadline and time.time() > deadline:
        return best['score'], best['state']

    key = normalize_key(state)

    # Base case or terminal: do level-0 resolution (random playout)
    if level == 0 or is_terminal(state, depth):
        if key in cache:
            score, s = cache[key]
        else:
            score, s = random_playout(state, all_moves, depth, ec_main)
            cache[key] = (score, s)
        if score > best['score']:
            best.update(score=score, state=s)
        return score, s

    # Enumerate legal moves (cap to e if provided)
    moves = legal_moves_fn(state, all_moves)
    if not moves:
        # No legal moves but not terminal: fall back to playout
        return _bilnmcs(state, 0, depth, all_moves, ec_main,ec_fast,
                      cache=cache, deadline=deadline, best=best,
                      tr=tr, trmax=trmax, b=b, r=r, e=e)

    if e is not None and len(moves) > e:
        moves = random.sample(moves, e)

    # Depth index (number of moves already chosen)
    d = min(len(state), depth)

    # ---------- (i) Cheap evaluation per move (b playouts) ----------
    candidates = []  # list of (mean_eval, move)
    for m in moves:
        if deadline and time.time() > deadline:
            return best['score'], best['state']

        s1 = normalize_sorted_list(list(state) + [m])
        k1 = normalize_key(s1)

        # Run b level-0 playouts from child to get a mean evaluation
        # We can cache per s1 playout result, but keep b>1 to reduce noise
        total = 0.0
        for _ in range(max(1, b)):
            # Use cache when luck helps, but still allow multiple draws
            if k1 in cache:
                sc, ss = cache[k1]
            else:
                sc, ss = random_playout(s1, all_moves, depth, ec_fast)
                cache[k1] = (sc, ss)
            total += sc
            # Track global best like NMCS
            if sc > best['score']:
                best.update(score=sc, state=ss)
        mean_eval = total / max(1, b)
        candidates.append((mean_eval, m))

        # ---------- (ii) Update per-depth running stats ----------
        # tr[d] keeps a running mean across all states at depth d
        acc = tr[d]
        acc['mean']  = (acc['mean'] * acc['count'] + mean_eval) / (acc['count'] + 1)
        acc['count'] += 1
        if mean_eval > trmax[d]:
            trmax[d] = mean_eval

    # ---------- (iii) Build depth-specific threshold and recurse ----------
    # θ_d = tr[d].mean + r * (trmax[d] - tr[d].mean)
    mu_d   = tr[d]['mean']
    best_d = trmax[d]
    theta  = mu_d + r * (best_d - mu_d)

    best_score_here, best_set_here = -float('inf'), []

    for mean_eval, m in candidates:
        if deadline and time.time() > deadline:
            return best['score'], best['state']

        s1 = normalize_sorted_list(list(state) + [m])
        k1 = normalize_key(s1)

        # Decide recursion depth: prune → level 0; else → level-1
        next_level = level - 1
        if mean_eval < theta:
            next_level = 0

        if k1 in cache and next_level == 0:
            score, s = cache[k1]
        else:
            score, s = _bilnmcs(s1, next_level, depth, all_moves, ec_main, ec_fast,
                              cache=cache, deadline=deadline, best=best,
                              tr=tr, trmax=trmax, b=b, r=r, e=e)
            # Only cache terminal/level-0 results, which are “complete” evaluations
            if next_level == 0:
                cache[k1] = (score, s)

        if score > best_score_here:
            best_score_here, best_set_here = score, s
            if score > best['score']:
                best.update(score=score, state=s)

    return best_score_here, best_set_here
'''

def bilnmcs(state, level, depth, all_moves, ec_main, ec_fast, *,
            b=2, r=0.5, e=None, timeout_sec=None):
    """
    Bi-Lazy NMCS entrypoint with per-level caches and expand-by-best loop.
    """
    deadline = time.time() + timeout_sec if timeout_sec else None

    # per-depthh running stats for thresholds
    max_depth = depth
    tr    = [{'mean': 0.0, 'count': 0} for _ in range(max_depth + 1)]
    trmax = [float('-inf') for _ in range(max_depth + 1)]

    # Per-level caches (main ans fast )
    c_main = [dict() for _ in range(level + 1)]
    c_fast = [dict() for _ in range(level + 1)]

    best = {'score': -float('inf'), 'state': []}
    score, s = _bilnmcs(state, level, depth, all_moves, ec_main, ec_fast,
                        c_main=c_main, c_fast=c_fast, deadline=deadline, best=best,
                        tr=tr, trmax=trmax, b=b, r=r, e=e)
    return score, normalize_sorted_list(s)


def _bilnmcs(state, level, depth, all_moves, ec_main, ec_fast, *,
             c_main, c_fast, deadline, best, tr, trmax, b, r, e):
    # Timeout 
    if deadline and time.time() > deadline:
        return best['score'], best['state']

    state = normalize_sorted_list(list(state))
    key   = normalize_key(state)

    if key in c_main[level]:
        return c_main[level][key]

    #  playout, cached at level 0
    if level == 0 or is_terminal(state, depth):
        if key in c_main[0]:
            score, s = c_main[0][key]
        else:
            score, s = random_playout(state, all_moves, depth, ec_main)
            c_main[0][key] = (score, s)
        if score > best['score']:
            best.update(score=score, state=s)
        return score, s

    bestSc_level = -float('inf')
    bestSet_level = []
    S_cur = list(state)  # we’ll grow this with next(bestSet \ S)

    while not is_terminal(S_cur, depth):
        if deadline and time.time() > deadline:
            return best['score'], best['state']

        # (i) Legal moves (cap to e)
        moves = legal_moves_fn(S_cur, all_moves)
        if not moves:
            break
        if e is not None and len(moves) > e:
            moves = random.sample(moves, e)

        # (ii) ccheap evals for each move using ec_fast, cached in c_fast[0]
        d = min(len(S_cur), depth)
        candidates = []  # (mean_eval, move)
        for m in moves:
            if deadline and time.time() > deadline:
                return best['score'], best['state']

            S1  = normalize_sorted_list(S_cur + [m])
            k1  = normalize_key(S1)
            tot = 0.0
            for _ in range(max(1, b)):
                if k1 in c_fast[0]:
                    sc, ss = c_fast[0][k1]
                else:
                    sc, ss = random_playout(S1, all_moves, depth, ec_fast)
                    c_fast[0][k1] = (sc, ss)
                tot += sc
                if sc > best['score']:
                    best.update(score=sc, state=ss)
            mean_eval = tot / max(1, b)
            candidates.append((mean_eval, m))

            # per-depth running stats
            acc = tr[d]
            acc['mean'] = (acc['mean'] * acc['count'] + mean_eval) / (acc['count'] + 1)
            acc['count'] += 1
            if mean_eval > trmax[d]:
                trmax[d] = mean_eval

        # Depth-specific threshold
        mu_d  = tr[d]['mean']
        b_d   = trmax[d]
        theta = mu_d + r * (b_d - mu_d)

        # Lazy recursion (prune -> level 0; else -> level-1), track best child-set
        local_best_sc, local_best_set = -float('inf'), []
        for mean_eval, m in candidates:
            if deadline and time.time() > deadline:
                return best['score'], best['state']

            S1 = normalize_sorted_list(S_cur + [m])
            k1 = normalize_key(S1)
            next_level = 0 if mean_eval < theta else (level - 1)

            if next_level == 0:
                if k1 in c_main[0]:
                    sc, s = c_main[0][k1]
                else:
                    sc, s = random_playout(S1, all_moves, depth, ec_main)
                    c_main[0][k1] = (sc, s)
            else:
                if k1 in c_main[level - 1]:
                    sc, s = c_main[level - 1][k1]
                else:
                    sc, s = _bilnmcs(S1, level - 1, depth, all_moves, ec_main, ec_fast,
                                     c_main=c_main, c_fast=c_fast, deadline=deadline, best=best,
                                     tr=tr, trmax=trmax, b=b, r=r, e=e)
                    c_main[level - 1][k1] = (sc, s)

            if sc > local_best_sc:
                local_best_sc, local_best_set = sc, s
            if sc > bestSc_level:
                bestSc_level, bestSet_level = sc, s
            if sc > best['score']:
                best.update(score=sc, state=s)

        # S_cur <- S_cur ∪ { next(bestSet \ S_cur) }
        x_star = next((x for x in local_best_set if x not in S_cur), None)
        if x_star is None:
            break  
        S_cur = normalize_sorted_list(S_cur + [x_star])
        
    return bestSc_level, S_cur
