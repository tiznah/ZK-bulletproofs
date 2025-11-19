from py_ecc.bn128 import G1, multiply, add, FQ, eq, Z1
from py_ecc.bn128 import curve_order as p
from functools import reduce
import random

def random_element():
    return random.randint(0, p - 1)

def add_points(*points):
    return reduce(add, points, Z1)

def vector_commit(points, scalars):
    if not points or not scalars: return Z1
    return reduce(add, [multiply(P, int(i)) for P, i in zip(points, scalars)], Z1)

G_vec = [
    (FQ(6286155310766333871795042970372566906087502116590250812133967451320632869759), FQ(2167390362195738854837661032213065766665495464946848931705307210578191331138)),
    (FQ(6981010364086016896956769942642952706715308592529989685498391604818592148727), FQ(8391728260743032188974275148610213338920590040698592463908691408719331517047)),
    (FQ(15884001095869889564203381122824453959747209506336645297496580404216889561240), FQ(14397810633193722880623034635043699457129665948506123809325193598213289127838)),
    (FQ(6756792584920245352684519836070422133746350830019496743562729072905353421352), FQ(3439606165356845334365677247963536173939840949797525638557303009070611741415))
]

H_vec = [
    (FQ(13728162449721098615672844430261112538072166300311022796820929618959450231493), FQ(12153831869428634344429877091952509453770659237731690203490954547715195222919)),
    (FQ(17471368056527239558513938898018115153923978020864896155502359766132274520000), FQ(4119036649831316606545646423655922855925839689145200049841234351186746829602)),
    (FQ(8730867317615040501447514540731627986093652356953339319572790273814347116534), FQ(14893717982647482203420298569283769907955720318948910457352917488298566832491)),
    (FQ(419294495583131907906527833396935901898733653748716080944177732964425683442), FQ(14467906227467164575975695599962977164932514254303603096093942297417329342836))
]

B = (FQ(11573005146564785208103371178835230411907837176583832948426162169859927052980), FQ(895714868375763218941449355207566659176623507506487912740163487331762446439))

# Folding functions
def fold_scalars(scalar_vec, u):
    u_inv = pow(u, p - 2, p)
    res = []
    for i in range(0, len(scalar_vec), 2):
        val = (scalar_vec[i] * u + scalar_vec[i+1] * u_inv) % p
        res.append(val)
    return res

def fold_points(point_vec, u):
    # Generators fold inversely: G_lo * u_inv + G_hi * u
    u_inv = pow(u, p - 2, p)
    res = []
    for i in range(0, len(point_vec), 2):
        pt = add(multiply(point_vec[i], u_inv), multiply(point_vec[i+1], u))
        res.append(pt)
    return res

def compute_LR(G, H, a, b, B_gen):

    # L terms
    L_terms_G = [multiply(G[i], a[i-1]) for i in range(1, len(G), 2)]
    L_terms_H = [multiply(H[i], b[i-1]) for i in range(1, len(H), 2)]
    
    # R terms
    R_terms_G = [multiply(G[i], a[i+1]) for i in range(0, len(G)-1, 2)]
    R_terms_H = [multiply(H[i], b[i+1]) for i in range(0, len(H)-1, 2)]
    
    # Blinding elements
    rho_L = random_element()
    rho_R = random_element()
    
    L = add_points(*L_terms_G, *L_terms_H, multiply(B_gen, rho_L))
    R = add_points(*R_terms_G, *R_terms_H, multiply(B_gen, rho_R))
    
    return L, R, rho_L, rho_R


def prove_logarithmic(a, b):
    n = len(a)
    
    G = G_vec[:n]
    H = H_vec[:n]
    alpha = random_element()
    
    # A = <a, G> + <b, H> + alpha * B
    comm_a = vector_commit(G, a)
    comm_b = vector_commit(H, b)
    comm_alpha = multiply(B, alpha)
    A = add_points(comm_a, comm_b, comm_alpha)
    

    curr_a, curr_b = a[:], b[:]
    curr_G, curr_H = G[:], H[:]
    curr_alpha = alpha
    curr_A = A
    
    round_num = 0
    while len(curr_a) > 1:
        if len(curr_a) % 2 != 0:
            print("Error: Vector length must be power of 2 for this implementation.")
            return False
            
        round_num += 1
        
        # Prover computes L, R (Cross-terms)
        L, R, rho_L, rho_R = compute_LR(curr_G, curr_H, curr_a, curr_b, B)

        # Interactive proof for simplicity        
        u = random_element() 
        
        # Folding
        curr_a = fold_scalars(curr_a, u)
        curr_b = fold_scalars(curr_b, u)
        
        # Update blinding
        u_sq = (u * u) % p
        u_inv_sq = pow(u_sq, p - 2, p)
        curr_alpha = (curr_alpha + u_sq * rho_L + u_inv_sq * rho_R) % p
        
        # Fold Generators
        curr_G = fold_points(curr_G, u)
        curr_H = fold_points(curr_H, u)
        
        # Update Commitment A for final check
        term_L = multiply(L, u_sq)
        term_R = multiply(R, u_inv_sq)
        curr_A = add_points(term_L, curr_A, term_R)
    
    final_a = curr_a[0]
    final_b = curr_b[0]
    final_alpha = curr_alpha
    final_G = curr_G[0]
    final_H = curr_H[0]
    
    LHS = curr_A 
    
    RHS_terms = [
        multiply(final_G, final_a),
        multiply(final_H, final_b),
        multiply(B, final_alpha)
    ]
    RHS = add_points(*RHS_terms)
    
    success = eq(LHS, RHS)
    assert success, "Proof Failed"
    print("Proof succeeded.")
    return

# Given tests
a1 = [808, 140, 166, 209]
b1 = [88, 242, 404, 602]

a2 = [433, 651]
b2 = [282, 521]

a3 = [222]
b3 = [313]

# Run proofs
prove_logarithmic(a1, b1)
prove_logarithmic(a2, b2)
prove_logarithmic(a3, b3)
