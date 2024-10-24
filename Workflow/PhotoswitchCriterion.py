import numpy as np

### script to analyze two spectra (open and close)
### for a given Spiropyran according to the
### Addressability criterion

#### helper functions
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

## gaussian applied on eV and not on nm 
def NmInEv(nm):
    eV_Val=1239.84193/nm
    return(eV_Val)


### main
def PhotoswitchCrit(SpecFormOne,SpecFormTwo,Routine='Quot',OpenorClosed='closed',NMCutoff=10,Wavelength_Cutoff=280,sigma=0.09):
    ### nm in eV everything will be converted in eV

    Wavelength_Cutoff=NmInEv(Wavelength_Cutoff)

    PhotoswitchCrit=0
    ### define UV and VIS form
    if OpenorClosed=='open':
        UVSpec=SpecFormTwo
        VisSpec=SpecFormOne
    if OpenorClosed=='closed':
        UVSpec=SpecFormOne
        VisSpec=SpecFormTwo

    UVSpec[:,1]=NmInEv(UVSpec[:,1])
    VisSpec[:,1]=NmInEv(VisSpec[:,1])

    ### define all closed states above 280 nm
    IndexList_UV=[]
    for i in range(len(UVSpec)):
        if UVSpec[i,1]<Wavelength_Cutoff:
            IndexList_UV.append(i)

    if len(IndexList_UV) == 0:
        return(0)


    ### Define x grid based on both spectra and define the full convoluted spectra
    X_Grid=list(UVSpec[:,1])+list(VisSpec[:,1])

    try:
        XRange=np.arange(np.round(np.min(X_Grid))+0.005,np.round(np.max(X_Grid))-0.005,0.001)
    except:
        return(PhotoswitchCrit)

    All_open_states=np.zeros(len(XRange))
    for i in range(len(VisSpec[:,1])):
        All_open_states=All_open_states+np.array(gaussian(XRange,VisSpec[i,1],sigma)*VisSpec[i,2])

    All_closed_states=np.zeros(len(XRange))
    for i in range(len(UVSpec[:,1])):
        All_closed_states=All_closed_states+np.array(gaussian(XRange,UVSpec[i,1],sigma)*UVSpec[i,2])

    if np.nan in XRange:
        return(0)

    if Routine=='Quot':
            ### define quot crit
            QuotCrit=np.zeros(len(XRange))
            for k in range(len(XRange)):
                if XRange[k] > Wavelength_Cutoff:
                    continue
                if All_closed_states[k]<0.05:
                    continue
                else:
                    QuotCrit[k]=All_closed_states[k]/(All_open_states[k]+0.05)

            return(np.max(QuotCrit))
