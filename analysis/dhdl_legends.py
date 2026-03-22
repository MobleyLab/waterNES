class Custom_Legends:
    def __init__(self):
        self.legend=[]
        self.legend1=[]
        self.legend2=[]
        self.legend3=[]
        self.legend4=[]

        self.legend.append(r'@ legend on')
        self.legend.append(r'@ legend box on')
        self.legend.append(r'@ legend loctype view')
        self.legend.append(r'@ legend 0.78, 0.8')
        self.legend.append(r'@ legend length 2')
        self.legend.append(r'@ s0 legend "Potential Energy (kJ/mol)"')
        self.legend.append(r'@ s1 legend "dH/d\xl\f{} vdw-lambda = 0.0000"')
        self.legend.append(r'@ s2 legend "dH/d\xl\f{} bonded-lambda = 0.0000"')
        self.legend.append(r'@ s3 legend "dH/d\xl\f{} restraint-lambda = 0.0000"')

        self.legend1=self.legend.copy()
        self.legend2=self.legend.copy()
        self.legend3=self.legend.copy()
        self.legend4=self.legend.copy()

        self.legend1.append(r'@ s4 legend "\xD\f{}H \xl\f{} to (0.0000, 0.0000, 0.0000)"')
        self.legend1.append(r'@ s5 legend "\xD\f{}H \xl\f{} to (0.0000, 1.0000, 0.0000)"')
        self.legend1.append(r'@ s6 legend "\xD\f{}H \xl\f{} to (1.0000, 1.0000, 0.0000)"')
        self.legend1.append(r'@ s7 legend "\xD\f{}H \xl\f{} to (1.0000, 1.0000, 0.0100)"')
        self.legend1.append(r'@ s8 legend "\xD\f{}H \xl\f{} to (1.0000, 1.0000, 0.0200)"')
        self.legend1.append(r'@ s9 legend "\xD\f{}H \xl\f{} to (1.0000, 1.0000, 0.0400)"')
        self.legend1.append(r'@ s10 legend "\xD\f{}H \xl\f{} to (1.0000, 1.0000, 0.0800)"')
        self.legend1.append(r'@ s11 legend "\xD\f{}H \xl\f{} to (1.0000, 1.0000, 0.1600)"')
        self.legend1.append(r'@ s12 legend "\xD\f{}H \xl\f{} to (1.0000, 1.0000, 0.3200)"')
        self.legend1.append(r'@ s13 legend "\xD\f{}H \xl\f{} to (1.0000, 1.0000, 0.6400)"')
        self.legend1.append(r'@ s14 legend "\xD\f{}H \xl\f{} to (1.0000, 1.0000, 1.0000)"')
        self.legend1.append(r'@ s15 legend "pV (kJ/mol)"')


        self.legend2.append(r'@ s4 legend "\xD\f{}H \xl\f{} to (0.0000, 0.0000, 0.0000)"')
        self.legend2.append(r'@ s5 legend "\xD\f{}H \xl\f{} to (0.0000, 1.0000, 0.0000)"')
        self.legend2.append(r'@ s6 legend "\xD\f{}H \xl\f{} to (1.0000, 1.0000, 0.0000)"')
        self.legend2.append(r'@ s7 legend "\xD\f{}H \xl\f{} to (1.0000, 1.0000, 1.0000)"')
        self.legend2.append(r'@ s8 legend "pV (kJ/mol)"')

        self.legend3.append(r'@ s4 legend "\xD\f{}H \xl\f{} to (0.0000, 0.0000, 0.0000)"')
        self.legend3.append(r'@ s5 legend "\xD\f{}H \xl\f{} to (1.0000, 1.0000, 1.0000)"')
        self.legend3.append(r'@ s6 legend "pV (kJ/mol)"')

        self.legend4.append(r'@ s4 legend "\xD\f{}H \xl\f{} to (0.0000, 0.0000, 0.0000)"')
        self.legend4.append(r'@ s5 legend "\xD\f{}H \xl\f{} to (0.0000, 1.0000, 0.0000)"')
        self.legend4.append(r'@ s6 legend "\xD\f{}H \xl\f{} to (1.0000, 1.0000, 1.0000)"')
        self.legend4.append(r'@ s7 legend "pV (kJ/mol)"')
"""
        if int(self.stage) < 5:
            if self.scheme == 1:
                return self.legend3
            elif self.scheme == 2:
                return self.legend2
            elif self.scheme == 3:
                return self.legend1
        elif int(self.stage) > 4:
            if self.scheme == 1:
                return self.legend2
            elif self.scheme == 2:
                return self.legend1

"""
