from dadi import Numerics, PhiManip, Integration, Spectrum

def simple_canary((nuIb0, nuIb, nuMo0, nuMo, nuCa0, TIbMo, TMoCa, TF),
    ns, pts):
    """
    nuIb0: starting Iberian population size
    nuIb: ending Iberian population size
    nuMo0: starting Moroccan population size
    nuMo: ending Morrocan population size
    nuCa0: starting Canary island population size
    nuCa: ending Canary island population size
    TIbMo: split time for Iberian and Moroccan populations
    TMoCa: split time for Moroccan and Canary island populations
    TF: simulation ending time
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, TIbMo, nu=nuIb0)
    print phi
    print TIbMo
    print xx
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, TMoCa, nu1=nuIb, nu2=nuMo0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, TF, nu1=nuIb, nu2=nuMo, nu3=nuCa0)


    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))

    return fs

def simple_canary_migration((nuIb0, nuIb, nuMo0, nuMo, nuCa0, TIbMo, TMoCa, mIbCa, TF),
    ns, pts):
    """
    nuIb0: starting Iberian population size
    nuIb: ending Iberian population size
    nuMo0: starting Moroccan population size
    nuMo: ending Morrocan population size
    nuCa0: starting Canary island population size
    nuCa: ending Canary island population size
    TIbMo: split time for Iberian and Moroccan populations
    TMoCa: split time for Moroccan and Canary island populations
    mIbCa: migration rate from Ib to Ca
    TF: simulation ending time
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, TIbMo, nu=nuIb0)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, TMoCa, nu1=nuIb, nu2=nuMo0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, TF, nu1=nuIb, nu2=nuMo, nu3=nuCa0, m13=mIbCa)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))

    return fs
