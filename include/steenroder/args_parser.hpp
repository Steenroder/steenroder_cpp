/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include <string>
#include <stdlib.h>
#include <stdexcept>

#include "anyoption.h"

namespace stn {

  class ArgsParser {
  private:
    AnyOption *option;

    void setUsage(AnyOption *opt) {
      opt->addUsage("Usage: ./steenroder [options] file...");
      opt->addUsage("Options: ");
      opt->addUsage(" -h  --help                     Display this information ");
      opt->addUsage(" -d  --dim <dim>                Dimension. Default: 1 ");
      opt->addUsage(" -k  --k <k>                    k. Default: 1 ");
      opt->addUsage(" -r  --reps                     Outputs representatives ");
      opt->addUsage("");
    }

    void setOptions(AnyOption *opt) {
      opt->setFlag("help", 'h');
      opt->setOption("dim", 'd');
      opt->setOption("k", 'k');
      opt->setFlag("reps", 'r');
    }

    AnyOption* initOption(int &argc, char **argv) {
      AnyOption *opt = new AnyOption();
      setUsage(opt);
      setOptions(opt);
      opt->processCommandArgs(argc, argv);
      return opt;
    }

    const char* getValue(char name, const char* default_value) {
      if (option->getValue(name) != NULL)
        return option->getValue(name);
      else
        return default_value;
    }

    const std::string getFilename(char* argv) {
      if (argv == NULL)
        throw std::runtime_error("No input file.");

      return argv;
    }


  public:
    const bool help;
    const unsigned int dim;
    const unsigned int k;
    const bool reps;
    const std::string input_filename;
    const std::string output_filename;

    ArgsParser(int &argc, char **argv)
      : option(initOption(argc, argv))
      , help(option->getFlag('h'))
      , dim(atoi(getValue('d', "1")))
      , k(atoi(getValue('k', "1")))
      , reps(option->getFlag('r'))
      , input_filename(getFilename(option->getArgv(0)))
      , output_filename(getFilename(option->getArgv(1)))
    {}

    ~ArgsParser() {
      delete option;
    }

    void printUsage() const {
        option->printUsage();
    }

  };

} // namespace stn
