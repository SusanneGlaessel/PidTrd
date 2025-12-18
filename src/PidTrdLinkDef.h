//#if defined(__CINT__) || defined(__CLING__)
//#ifdef __MAKECINT__
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class TrdContainer +;
#pragma link C++ class PidTrd::ParticleProb + ;
#pragma link C++ class PidTrd::Getter + ;


#endif
