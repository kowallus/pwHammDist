#include "helper.h"

#include "byteswap.h"

std::ostream *logout = &std::cout;

void cleanCache() {
    const int ccSize = 10*1024*1024; // Allocate 10M. Set much larger then L2
    const int ccRepeats = 0xff;
    char *cc = (char *)malloc(ccSize);
    for (int cci = 0; cci < ccRepeats; cci++)
        for (int ccj = 0; ccj < ccSize; ccj++)
            cc[ccj] = cci+ccj;
    free(cc);
}

// TIME

clock_t checkpoint;

clock_t clock_checkpoint() {
    checkpoint = clock();
    return checkpoint;
//    cout << "Clock reset!\n";
}

unsigned long long int clock_millis(clock_t checkpoint) {
    return (clock() - checkpoint) * (unsigned long long int) 1000 / CLOCKS_PER_SEC;
}

unsigned long long int clock_millis() {
    return clock_millis(checkpoint);
}

chrono::steady_clock::time_point chronocheckpoint;

chrono::steady_clock::time_point time_checkpoint() {
    chronocheckpoint = chrono::steady_clock::now();
    return chronocheckpoint;
}

unsigned long long int time_millis(chrono::steady_clock::time_point checkpoint) {
    chrono::nanoseconds time_span = chrono::duration_cast<chrono::nanoseconds>(chrono::steady_clock::now() - checkpoint);
    return (double)time_span.count() / 1000000.0;
}

unsigned long long int time_millis() {
    return time_millis(chronocheckpoint);
}

const size_t chunkSize = 10000000;

void* readArray(std::istream& in, size_t arraySizeInBytes) {

    if (in) {
        char* destArray = new char[arraySizeInBytes];
        readArray(in, destArray, arraySizeInBytes);
        return (void*) destArray;
    } else
        throw(errno);
}

void readArray(std::istream& in, void* destArray, size_t arraySizeInBytes) {

    if (in) {
        size_t length = arraySizeInBytes;
        char* ptr = (char*) destArray;
        size_t bytesLeft = length;
        while (bytesLeft > chunkSize) {
            in.read(ptr, chunkSize);
            ptr = ptr + chunkSize;
            bytesLeft -= chunkSize;
        }
        in.read(ptr, bytesLeft);
    } else
        throw(errno);
}

void writeArray(std::ostream& out, void* srcArray, size_t arraySize, bool verbose) {
    size_t bytesLeft = arraySize;
    if (out) {
        while (bytesLeft > chunkSize) {
            out.write((char*) srcArray, chunkSize);
            srcArray = (void*) (((char*) srcArray) + chunkSize);
            bytesLeft -= chunkSize;
        }
        out.write((char*) srcArray, bytesLeft);

    } else
        throw(errno);
    if (verbose)
        cout << "Written " << arraySize << " bytes\n";
}

void* readWholeArray(std::istream& in, size_t& arraySizeInBytes) {

    if (in) {
        in.seekg(0, in.end);
        size_t length = in.tellg();
        char* destArray = new char[length];
        in.seekg(0, in.beg);
        char* ptr = destArray;
        size_t bytesLeft = length;
        while (bytesLeft > chunkSize) {
            in.read(ptr, chunkSize);
            ptr = ptr + chunkSize;
            bytesLeft -= chunkSize;
        }
        in.read(ptr, bytesLeft);

        arraySizeInBytes = length;
        return (void*) destArray;
    } else
        throw(errno);
}

void* readWholeArrayFromFile(string srcFile, size_t& arraySizeInBytes) {
    clock_checkpoint();

    std::ifstream in(srcFile.c_str(), std::ifstream::binary);

    void* destArray = readWholeArray(in, arraySizeInBytes);

    cout << "Read " << arraySizeInBytes << " bytes from " << srcFile << " in " << clock_millis() << " msec \n";

    return destArray;
}

void writeArrayToFile(string destFile, void* srcArray, size_t arraySize) {
    clock_checkpoint();

    std::ofstream out(destFile.c_str(), std::ios::out | std::ios::binary);

    writeArray(out, srcArray, arraySize);

    cout << "Write " << arraySize << " bytes to " << destFile << " in " << clock_millis() << " msec \n";
}

void writeStringToFile(string destFile, const string &src) {
    writeArrayToFile(destFile, (void*) src.data(), src.length());
}

bool plainTextWriteMode = false;

void writeReadMode(std::ostream &dest, bool plainTextWriteMode) {
    dest << (plainTextWriteMode?TEXT_MODE_ID:BINARY_MODE_ID) << "\n";
}

bool confirmTextReadMode(std::istream &src) {
    string readMode;
    src >> readMode;
    if (readMode != TEXT_MODE_ID && readMode != BINARY_MODE_ID) {
        fprintf(stderr, "Expected READ MODE id (not: %s)\n", readMode.c_str());
        exit(EXIT_FAILURE);
    }
    char check = src.get();
    return readMode == TEXT_MODE_ID;
}

template<>
void writeValue(std::ostream &dest, const uint8_t value, bool plainTextWriteMode) {
    if (plainTextWriteMode)
        dest << (uint16_t) value << endl;
    else
        dest.write((char *) &value, sizeof(uint8_t));
}

template<>
void readValue(std::istream &src, uint8_t &value, bool plainTextReadMode) {
    if (plainTextReadMode) {
        uint16_t temp;
        src >> temp;
        value = (uint8_t) temp;
    } else
        src.read((char *) &value, sizeof(uint8_t));
}

void writeUIntByteFrugal(std::ostream &dest, uint64_t value) {
    while (value >= 128) {
        uint8_t yByte = 128 + (value % 128);
        dest.write((char *) &yByte, sizeof(uint8_t));
        value = value / 128;
    }
    uint8_t yByte = value;
    dest.write((char *) &yByte, sizeof(uint8_t));
}

bool bytePerReadLengthMode = false;

void readReadLengthValue(std::istream &src, uint16_t &value, bool plainTextReadMode) {
    if (bytePerReadLengthMode)
        readValue<uint8_t>(src, (uint8_t&) value, plainTextReadMode);
    else
        readValue<uint16_t>(src, value, plainTextReadMode);
}

void writeReadLengthValue(std::ostream &dest, const uint16_t value) {
    if (bytePerReadLengthMode)
        writeValue<uint8_t>(dest, (uint8_t) value, plainTextWriteMode);
    else
        writeValue<uint16_t>(dest, value, plainTextWriteMode);
}

string toString(unsigned long long value) {
        std::ostringstream oss;
        oss << value;
        return oss.str();
};

string toMB(unsigned long long value, unsigned char decimalPlaces) {
    std::ostringstream oss;
    int power = 1000000;
    oss << value / power;
    if (decimalPlaces > 0) {
        oss << ".";
        for(int i = 0; i < decimalPlaces; i++) 
            oss << ((value / (power /= 10)) % 10);
        
    }
    return oss.str();
}

string toString(long double value, unsigned char decimalPlaces) {
    std::ostringstream oss;
    oss << (long long int) value;
    if (decimalPlaces > 0) {
        oss << ".";
        int power = 1000000;        
        unsigned long long decimals = (value - (long long int) value) * power;
        for(int i = 0; i < decimalPlaces; i++) 
            oss << ((decimals / (power /= 10)) % 10);
    }
    return oss.str();
}

string alignRight(string str, unsigned char width) {
    int spaces = width - str.length();
    while(spaces-- > 0)
        str.insert(str.begin(), ' ');
    return str;
};

unsigned long long int powuint(unsigned long long int base, int exp)
{
    if (exp == 0) return 1;
    if (exp == 1) return base;

    unsigned long long int tmp = powuint(base, exp/2);
    if (exp%2 == 0) return tmp * tmp;
        else return base * tmp * tmp;
}

int readsSufPreCmp(const char* suffixPart, const char* prefixRead) {
    while (*suffixPart) {
        if (*suffixPart > *prefixRead)
            return 1;
        if (*suffixPart++ < *prefixRead++)
            return -1;
    }
    return 0;
}

int strcmplcp(const char* lStrPtr, const char* rStrPtr, int length) {
    
    int i = 0;
    while (length - i >= 4) {
        int cmp = bswap_32(*(uint32_t*) lStrPtr) - bswap_32(*(uint32_t*) rStrPtr);
        if (cmp != 0)
            break;
        lStrPtr += 4;
        rStrPtr += 4;
        i += 4;
    }

    while (i < length) {
        i++;
        int cmp = *(unsigned char*)(lStrPtr++) - *(unsigned char*)(rStrPtr++);
        if (cmp > 0)
            return 1;
        if (cmp < 0)
            return -1;
    }

    return 0;

}
