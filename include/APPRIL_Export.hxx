#ifndef ARP_EXPORT_HXX_INCLUDED
#define ARP_EXPORT_HXX_INCLUDED

#if defined _WINDOWS || WINDOWS

	#ifdef APPRIL_VE_EXPORTS

	#define APPRIL_EXP_FUNC  __declspec(dllexport)
	#define APPRIL_EXP_CLASS  __declspec(dllexport)
	#define APPRIL_EXP_OBJ  __declspec(dllexport)
	#define APPRIL_TEMPLATE_INSTANTIATION

	#else

	#define APPRIL_EXP_FUNC  __declspec(dllimport)
	#define APPRIL_EXP_CLASS  __declspec(dllimport)
	#define APPRIL_EXP_OBJ  __declspec(dllimport)
	#define APPRIL_TEMPLATE_INSTANTIATION extern

	#endif

#else

	#define APPRIL_EXP_FUNC  
	#define APPRIL_EXP_CLASS 
	#define APPRIL_EXP_OBJ  
	#define APPRIL_TEMPLATE_INSTANTIATION

#endif

#endif // ARP_EXPORT_HXX_INCLUDED
