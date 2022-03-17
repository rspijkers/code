# particle class. 
# Example: Kplus = hadron('K+', 321, "K^{+}")
class hadron:
    def __init__(self, name: str, pdg: int, latex: str) -> None:
        # note that 
        self.name = name
        self.pdg = pdg
        self.latex = latex

# define here the particles, try to group them in a pseudological fashion:
# particles start with a capital letter; greek letters and +/0/- are written out
# for now only strange hadrons are written out
Kzero = hadron('K0', 311, "K^{0}")
Kzerobar = hadron('K0bar', -311, "#bar{K^{0}}")
Kplus = hadron('K+', 321, 'K^{+}')
Kminus = hadron('K-', -321, 'K^{-}')

Lambda = hadron('Lambda', 3122, '#Lambda') 
Lambdabar = hadron('Lambdabar', -3122, '#bar{#Lambda}')

Sigmazero = hadron('Sigma0', 3212, '#Sigma^{0}')
Sigmazerobar = hadron('Sigma0bar', -3212, '#bar{#Sigma^{0}}')
Sigmaminus = hadron('Sigma-', 3112, '#Sigma^{-}')
Sigmaminusbar = hadron('Sigma-bar', -3112, '#bar{#Sigma^{-}}')
Sigmaplus = hadron('Sigma+', 3222, '#Sigma^{+}')
Sigmaplusbar = hadron('Sigma+bar', -3222, '#bar{#Sigma^{+}}')

Xizero = hadron('Xi0', 3322, '#Xi^{0}')
Xizerobar = hadron('Xi0bar', -3322, '#bar{#Xi^{0}}')
Ximinus = hadron('Xi-', 3312, '#Xi^{-}')
Xiplus = hadron('Xi+', -3312, '#Xi^{+}')

Omegaminus = hadron('Omega-', 3334, '#Omega^{-}')
Omegaplus = hadron('Omega+', -3334, '#Omega^{+}')
