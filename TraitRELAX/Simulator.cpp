#include <TraitRELAXManager.h>


int main(int args, char** argv)
{
  cout << "************************************************************************************************" << endl;
  cout << "*                                     TraitRELAX simulator                                     *" << endl;
  cout << "************************************************************************************************" << endl;
  cout << endl;

  try
  { 
    // process input from params file
    BppApplication simulationParams(args, argv, "simulationParams");

    TraitRELAXManager* traitRELAXManager = new TraitRELAXManager(&simulationParams);

    traitRELAXManager->simulate();

    // free
    delete traitRELAXManager;
  }
  
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}