import math

def fourvec(pt, eta, phi, mass=0.0):
    px = pt * math.cos(phi)
    py = pt * math.sin(phi)
    pz = pt * math.sinh(eta)
    E  = math.sqrt(px*px + py*py + pz*pz + mass*mass)
    return (E, px, py, pz)

def add_fv(a, b):
    return tuple(x+y for x,y in zip(a,b))

def mass_from_fv(fv):
    E, px, py, pz = fv
    m2 = E*E - (px*px + py*py + pz*pz)
    return math.sqrt(m2) if m2>0 else 0.0

def eta_phi_from_fv(fv):
    E, px, py, pz = fv
    pt = math.hypot(px, py)
    phi = math.atan2(py, px)
    p = math.sqrt(px*px + py*py + pz*pz)
    # protect against division by zero
    if p == abs(pz):
        eta = float("inf") if pz>0 else float("-inf")
    else:
        eta = 0.5 * math.log((p + pz) / (p - pz))
    return eta, phi, pt

def check_pair(parent, g1, g2):
    # parent, g1, g2 are dict-like with keys: pt, eta, phi, mass (mass optional)
    fv_p = fourvec(parent["pt"], parent["eta"], parent["phi"], parent.get("mass", 0.0))
    fv1  = fourvec(g1["pt"], g1["eta"], g1["phi"])
    fv2  = fourvec(g2["pt"], g2["eta"], g2["phi"])
    fv_sum = add_fv(fv1, fv2)

    m_gg = mass_from_fv(fv_sum)
    eta_sum, phi_sum, pt_sum = eta_phi_from_fv(fv_sum)

    # deltaR between sum direction and parent
    dphi = phi_sum - parent["phi"]
    while dphi > math.pi: dphi -= 2*math.pi
    while dphi < -math.pi: dphi += 2*math.pi
    dR = math.sqrt((eta_sum - parent["eta"])**2 + dphi**2)

    rel_pt = (pt_sum - parent["pt"]) / parent["pt"]
    rel_m  = (m_gg - parent.get("mass", m_gg)) / parent.get("mass", m_gg)

    return {
        "m_gg": m_gg,
        "pt_sum": pt_sum,
        "deltaR": dR,
        "rel_pt_diff": rel_pt,
        "rel_mass_diff": rel_m
    }

# Example (use your numbers)
parent1 = {"pt":19.38, "eta":-3.89, "phi":-0.85, "mass":15.0}
g1 = {"pt":15.62,  "eta":-3.90, "phi":-1.23}
g2 = {"pt":4.75, "eta":0.09, "phi":1.0}
print(check_pair(parent1, g1, g2))

parent2 = {"pt":7.12, "eta":-0.07, "phi":-1.91, "mass":15.00}
g3 = {"pt":7.55, "eta":-3.18, "phi":0.02}
g4 = {"pt":11.81, "eta":-0.08, "phi":-2.00}
print(check_pair(parent2, g3, g4))
