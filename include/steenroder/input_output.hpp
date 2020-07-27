/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include <sys/stat.h>
#include <iostream>
#include <string>

#include "commons.hpp"
#include "options.hpp"

namespace stn {

  template <class T, InputOutputFormat inputOutput, InputOutputMode inputOutputMode>
    class InputOutput {};

  template <class T>
    class InputOutput<T, InputOutputFormat::Generic, InputOutputMode::Generic> {
  protected:
    const std::string writerFolder;
    const std::string fileExtension;
    const std::string filePrefix;
    const std::string fileMode;

    std::ofstream file;

    InputOutput(const std::string& writerFolder_in,
           const std::string& filePrefix_in,
           const std::string& fileExtension_in,
           const std::string& fileMode_in)
      : writerFolder(writerFolder_in)
      , fileExtension(fileExtension_in)
      , filePrefix(filePrefix_in)
      , fileMode(fileMode_in)
    {}

    inline std::string getFileName(const std::string& postfix = "") {
      return writerFolder + filePrefix + postfix + fileExtension;
    }

    inline void closeFile() { Base::file.close(); }

  };


  template <class T>
    class InputOutput<T, InputOutputFormat::Generic, InputOutputMode::ascii>
    : public InputOutput<T, InputOutputFormat::Generic, InputOutputMode::Generic> {
  private:
    using Base = InputOutput<T, InputOutputFormat::Generic, InputOutputMode::Generic>;

  public:
    InputOutput(const std::string& writerFolder_in,
                const std::string& filePrefix_in,
                const std::string& fileExtension_in)
      : Base(writerFolder_in + "/", filePrefix_in, fileExtension_in, "ascii") {}

    using Base::getIsWritten;

  protected:
    inline void openAndAppend(const std::string& fileName) {
      Base::file.open(fileName, std::ofstream::out | std::ofstream::app);
      Base::file.precision(16);

      if (!Base::file) {
        std::cout << "Could not open file " << fileName << std::endl;
      }
    }

    inline void openAndTruncate(const std::string& fileName) {
      Base::file.open(fileName, std::ofstream::out | std::ofstream::trunc);
      Base::file.precision(16);

      if (!Base::file) {
        std::cout << "Could not open file " << fileName << std::endl;
      }
    }

    template <class U>
    inline void write(const U data) {
      Base::file << data;
    }

  };


  template <class T>
  class InputOutput<T, InputOutputFormat::Generic, InputOutputMode::binary>
    : public InputOutput<T, InputOutputFormat::Generic, InputOutputMode::Generic> {
  private:
    using Base = InputOutput<T, InputOutputFormat::Generic, InputOutputMode::Generic>;

  public:
    InputOutput(const std::string& writerFolder_in,
                const std::string& filePrefix_in,
                const std::string& fileExtension_in)
      : Base(writerFolder_in + "/", filePrefix_in, fileExtension_in, "binary")
    {}

    using Base::getIsWritten;

  protected:
    inline void openAndAppend(const std::string& fileName) {
      Base::file.open(fileName, std::ofstream::out | std::ofstream::app |
                      std::ofstream::binary);
    }

    inline void openAndTruncate(const std::string& fileName) {
      Base::file.open(fileName, std::ofstream::out | std::ofstream::trunc |
                      std::ofstream::binary);
    }

    template <class U>
    inline void write(U data) {
      Base::file.write(reinterpret_cast<char*>(&data), sizeof(data));
    }
  };

}  // namespace stn
