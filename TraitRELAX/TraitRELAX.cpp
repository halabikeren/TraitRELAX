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

    if(args == 1)
    {
      cout << "Please provide parameter file. Refer to https://github.com/halabikeren/TraitRELAX/Examples/example2.bpp for example" << endl;
      return 0;
    }

    BppApplication* traitRELAXParameters = new BppApplication(args, argv, "traitRELAX");

    TraitRELAXManager* traitRELAXManager = new TraitRELAXManager(traitRELAXParameters);
    
    traitRELAXManager->init();

    map<string, double> optimizedNullModelParameters = traitRELAXManager->optimizeNullModel();

    map<string, double> optimizedAlternativeModelParameters = traitRELAXManager->optimizeAlternativeModel();

    traitRELAXManager->test(optimizedNullModelParameters, optimizedAlternativeModelParameters); // Perform a statistical test

    delete traitRELAXManager;
    delete traitRELAXParameters;
  }
  catch (exception &e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}