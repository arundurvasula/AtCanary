from dadi import Numerics, PhiManip, Integration, Spectrum

def simple_canary((nuIb0, nuIb, nuMo0, nuMo, nuCa0, TIbMo, TMoCa, TF),
    ns, pts):
    """
    nuIb0: starting Iberian population size
    nuIb: ending Iberian population size
    nuMo0: starting Moroccan population size
    nuMo: ending Morrocan population size
    nuCa0: starting Canary island population size
    TIbMo: split time for Iberian and Moroccan populations
    TMoCa: split time for Moroccan and Canary island populations
    TF: simulation ending time
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, TIbMo, nu=nuIb0)

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

def two_step((nMo, nCa, T1, T2, T3), ns, pts):
    """
    nMo: Moroccan constant population size
    nCa: Canary Island constant population size
    T1: time for first split (Iberia to Iberia and Morocco)
    T2: time for second split (Morocco to Morocco and Canary)
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, T1, nu=1)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T2, nu1=1, nu2=nMo)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T3, nu1=1, nu2=nMo, nu3=nCa)

    fs = Spectrum.from_phi(phi, ns, (xx, xx, xx))

    return fs

def two_step_migration((nMo, nCa, mIbCa, T1, T2), ns, pts):
    """
    nMo: Moroccan constant population size
    nCa: Canary Island constant population size
    T1: time for first split (Iberia to Iberia and Morocco)
    T2: time for second split (Morocco to Morocco and Canary)
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1=1, nu2=nMo)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T2, nu1=1, nu2=nMo, nu3=nCa, m13=mIbCa)

    fs = Spectrum.from_phi(phi, ns, (xx, xx, xx))

    return fs

def Mo_Ib((nIb, nMo1, nMo2, T1), ns, pts):
    """
    nIb: Iberian constant population size
    nMo1: starting Moroccan size
    nMo2: ending Moroccan size
    T1: Time to first split
    This model doesn't consider the Canary Island population size or split time.
    """

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1=nIb, nu2=nMo1)

    nu_funcMo = lambda t: (nMo2-1) * (t/T1) + 1

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, 1, nu1=nIb, nu2=nu_funcMo, nu3=1, m13=0.2)

    fs = Spectrum.from_phi(phi, ns, (xx, xx, xx))

    return fs

def moroccoIberia_migration_Andrea((nuIb, nuMo, TIbMo, mIbMo, mMoIb), 
    ns, pts):
    """
    nuIb: Iberian population size
    nuMo: Moroccan population size
    TIbMo: split time for Iberian and Moroccan populations
    mIbmo: migration rate from Ib to Mo
    mMoIb: migration rate from Mo to Ib
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, TIbMo, nu1=nuIb, nu2=nuMo, m12=mIbMo, m21=mMoIb)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))

    return fs

def moroccoIberia_migration((nuIb, nuMo, nuCa, TIbMo, TMoCa, mIbMo, mMoIb, mMoCa),
    ns, pts):
    """
    nuIb: Iberian population size
    nuMo: Moroccan population size
    nuCa: Canary population size
    TIbMo: split time for Iberian and Moroccan populations
    TMoCa: split time for Moroccan and Canary populations
    mIbmo: migration rate from Ib to Mo
    mMoIb: migration rate from Mo to Ib
    mMoCa: migration rate from Mo to Ca
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, TIbMo, nu1=nuIb, nu2=nuMo, m12=mIbMo, m21=mMoIb)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, TMoCa, nu1=nuIb, nu2=nuMo, nu3=nuCa, m13=mMoCa)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))

    return fs
