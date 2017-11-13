#include "Problems/ch1_cantilever.hh"
#include "Problems/ch1_helixTerminalMoment.hh"
#include "Problems/ch5_helixUncoiling.hh"
#include "Problems/ch8_vibratingCable.hh"

void CreateProblemVector()
{
  problems.push_back(NULL);
  problems.push_back(new ch1_cantilever());
  problems.push_back(new ch1_helixTerminalMoment());
  problems.push_back(new ch5_helixUncoiling());
  problems.push_back(new ch8_vibratingCable());
}
