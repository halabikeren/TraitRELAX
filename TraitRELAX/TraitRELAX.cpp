#include <TraitRELAXManager.h>

/******************************************************************************/
/*********************************** Main *************************************/
/******************************************************************************/

int main(int args, char **argv)
{
  cout << "************************************************************************************************" << endl;
  cout << "*       TraitRELAX program for detecting phenotype related changes in selective pressure       *" << endl;
  cout << "*                                                                                              *" << endl;
  cout << "*                       by Halabi, Levy Karin, Guéguen and Mayrose 2020                        *" << endl;
  cout << "************************************************************************************************" << endl;
  cout << endl;

  try
  {
    BppApplication traitRELAXParameters(args, argv, "traitRELAX");

    TraitRELAXManager* traitRELAXManager = new TraitRELAXManager(&traitRELAXParameters);
    
    traitRELAXManager->init();

    map<string, double> optimizedNullModelParameters = traitRELAXManager->optimizeNullModel();

    map<string, double> optimizedAlternativeModelParameters = traitRELAXManager->optimizeAlternativeModel();

    traitRELAXManager->test(optimizedNullModelParameters, optimizedAlternativeModelParameters); // Perform a statistical test

    delete traitRELAXManager;
  }
  catch (exception &e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}