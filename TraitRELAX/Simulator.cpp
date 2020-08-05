#include <TraitRELAXManager.h>


int main(int args, char** argv)
{
  cout << "************************************************************************************************" << endl;
  cout << "*                                     TraitRELAX simulator                                     *" << endl;
  cout << "************************************************************************************************" << endl;
  cout << endl;

  try
  { 

    if(args == 1)
    {
      cout << "Please provide parameter file. Refer to https://github.com/halabikeren/TraitRELAX/Examples/simulate.bpp for example" << endl;
      return 0;
    }

    // process input from params file
    BppApplication* simulationParams = new BppApplication(args, argv, "simulationParams");

    TraitRELAXManager* traitRELAXManager = new TraitRELAXManager(simulationParams);

    traitRELAXManager->simulate();

    // free
    delete traitRELAXManager;
    delete simulationParams;
  }
  
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}