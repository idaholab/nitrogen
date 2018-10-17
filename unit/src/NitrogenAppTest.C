#include <gtest/gtest.h>
#include <dlfcn.h>
#include "pcrecpp.h"
#include "NitrogenApp.h"
#include "Factory.h"
#include "AppFactory.h"
#include "Syntax.h"

extern char ** my_argv;

static const std::string AppName = "NitrogenApp";
static const std::string LibName = "libnitrogen";
static const std::string LibPath = "../lib";

/**
 * Get the library name from the binary we are running
 * @param[in] binary_filename The name of the binary we are running
 * @return The file name of the library that contains the entry points we will need later
 */
std::string
getLibName(const std::string & binary_filename)
{
  pcrecpp::RE re_method(".+-(.+)$");
  std::string method;
  if (re_method.FullMatch(binary_filename, &method))
    return LibName + "-" + method + ".la";
  else
    return "";
}

/**
 * Extract dlname from the .la file
 * .la is portable and will have different lib name on loinux and MacOS X
 * @param library_filename The name of the .la file
 * @return The content of dlname found in the .la file
 */
std::string
getDlname(const std::string & library_filename)
{
  std::string line;
  std::string dl_lib_filename;

  std::ifstream handle(library_filename.c_str());
  if (handle.is_open())
  {
    while (std::getline(handle, line))
    {
      // Look for the system dependent dynamic library filename to open
      if (line.find("dlname=") != std::string::npos)
        // Magic numbers are computed from length of this string "dlname=' and line minus that
        // string plus quotes"
        dl_lib_filename = line.substr(8, line.size() - 9);
    }
    handle.close();
  }

  return dl_lib_filename;
}

TEST(NitrogenApp, external_registerApp)
{
  std::string libname = getLibName(my_argv[0]);
  ASSERT_TRUE(libname.size() > 0);

  std::string dlname = LibPath + "/" + getDlname(LibPath + "/" + libname);
  ASSERT_TRUE(dlname.size() > 0);

#ifdef LIBMESH_HAVE_DLOPEN
  void * handle = dlopen(dlname.c_str(), RTLD_LAZY);
  void * registration_method = dlsym(handle, std::string(AppName + "__registerApps").c_str());

  if (registration_method)
  {
    typedef void (*register_app_t)();
    register_app_t * reg_ptr = reinterpret_cast<register_app_t *>(&registration_method);
    (*reg_ptr)();

    int argc = 0;
    char * argv[1] = {NULL};
    std::shared_ptr<MooseApp> app = AppFactory::createAppShared(AppName, argc, argv);
    ASSERT_TRUE(app.get() != nullptr);
  }
  else
    ASSERT_TRUE(false);
  dlclose(handle);
#endif
}

TEST(NitrogenApp, external_registerAll)
{
  std::string libname = getLibName(my_argv[0]);
  ASSERT_TRUE(libname.size() > 0);

  std::string dlname = LibPath + "/" + getDlname(LibPath + "/" + libname);
  ASSERT_TRUE(dlname.size() > 0);

#ifdef LIBMESH_HAVE_DLOPEN
  void * handle = dlopen(dlname.c_str(), RTLD_LAZY);
  void * registration_method = dlsym(handle, std::string(AppName + "__registerAll").c_str());

  if (registration_method)
  {
    typedef void (*register_app_t)(Factory *, ActionFactory * af, Syntax * s);
    register_app_t * reg_ptr = reinterpret_cast<register_app_t *>(&registration_method);

    std::shared_ptr<MooseApp> app = AppFactory::createAppShared(AppName, 0, nullptr);
    (*reg_ptr)(&app->getFactory(), &app->getActionFactory(), &app->syntax());
  }
  else
    ASSERT_TRUE(false);
  dlclose(handle);
#endif
}
