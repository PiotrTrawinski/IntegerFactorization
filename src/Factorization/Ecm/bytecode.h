#pragma once
#include <cstdint>
#include <vector>
#include <map>

namespace bytecode {
    enum class Block {
        End = 0,
        Naf = 1,
        DbChain = 2,
        Prac = 11,
    };
    
    enum class NafOpCode : uint8_t {
        END  = 0b000,
        ADDn = 0b001,
        SUBn = 0b010,
        DBLn = 0b011,

        FullMask = 0b100,
        ADD  = 0b101,
        SUB  = 0b110,
        DBL  = 0b111,
    };

    enum class PracOpCode : uint8_t {
        Rule1 = 1,
        Rule2 = 2,
        Rule3 = 3,
        Rule4 = 4,
        Rule5 = 5,
        Rule6 = 6,
        Rule7 = 7,
        Rule8 = 8,
        Rule9 = 9,
        End = 10,
    };

    /*
        File Format:
        64 bits: size in bytes
        64 bits: B
        32 bits: DBL  count
        32 bits: DBLn count
        32 bits: TPL  count
        32 bits: TPLn count
        32 bits: ADD  count
        32 bits: ADDn count
        32 bits: Dadd count
        32 bits: DDbl count
        ?? bits: bytecode
    */

    struct Buffer {
        std::vector<uint8_t> data;
        int size = 0;

        Buffer(int capacity=100'000) {
            data.resize(capacity);
        }
        void push(uint8_t n) {
            data[size++] = n;
        }
        void push32(uint32_t n) {
            push(n);
            push(n >> 8);
            push(n >> 16);
            push(n >> 24);
        }
        void push64(uint64_t n) {
            push32(n);
            push32(n >> 32);
        }
        void insert32(int index, uint32_t n) {
            data[index]     = n;
            data[index + 1] = n >> 8;
            data[index + 2] = n >> 16;
            data[index + 3] = n >> 24;
        }
        void enforceCapacityBigEnough() {
            if (size + 50'000 >= data.capacity()) {
                data.resize(data.capacity() * 2);
            }
        }
    };

    struct Instruction {
        uint8_t dblCount;
        uint8_t tplCount;
        uint8_t index;
        bool isSub;
        bool skipAdd;
        bool isFinal;
    };

    struct Reader {
        const std::vector<uint8_t>& buffer;
        int readIndex = 0;

        Reader(const std::vector<uint8_t>& buffer, int startPos=0) : buffer(buffer) {
            readIndex = startPos;
        }

        void skipByte() {
            readIndex++;
        }
        Block peekNextBlockOpCode() {
            return static_cast<Block>(buffer[readIndex] & 0x0f);
        }
        NafOpCode peekNafOpCode() {
            return static_cast<NafOpCode>(buffer[readIndex] & 0x0f);
        }
        bool peekIfPracSwap() {
            return buffer[readIndex] & 0x10;
        }
        int peekRepCount() {
            return ((buffer[readIndex] & 0b11100000) >> 5) + 1;
        }
        PracOpCode peekPracOpCode() {
            return static_cast<PracOpCode>(buffer[readIndex] & 0x0f);
        }
        uint8_t peekDataBits() {
            return buffer[readIndex] >> 4;
        }

        int getPointIndex() {
            return buffer[readIndex++] & 0x0f;
        }
        int getNextByte() {
            return buffer[readIndex++];
        }
        uint8_t peekByte() {
            return buffer[readIndex];
        }

        Instruction nextInstruction() {
            Instruction inst;
            inst.skipAdd  = buffer[readIndex] & 0x40;
            inst.isFinal  = buffer[readIndex] & 0x20;
            inst.isSub    = buffer[readIndex] & 0x10;
            inst.index    = buffer[readIndex] & 0x0f;
            inst.dblCount = buffer[readIndex + 1];
            if (buffer[readIndex] >> 7) {
                inst.tplCount = buffer[readIndex + 2];
                readIndex += 3;
            } else {
                inst.tplCount = 0;
                readIndex += 2;
            }
            return inst;
        }
    };

    struct OperationCounts {
        uint32_t dbl  = 0;
        uint32_t dbln = 0;
        uint32_t tpl  = 0;
        uint32_t tpln = 0;
        uint32_t add  = 0;
        uint32_t addn = 0;
        uint32_t dadd = 0;
        uint32_t ddbl = 0;
    };

    struct Writer {
        Buffer buffer;
        uint8_t pracLastByte = 0;
        int curDblCount = 0;
        int curTplCount = 0;
        int lastInstructionId = 0;
        int sizePosition = 0;
        OperationCounts operCounts;

        void START(uint64_t B) {
            pracLastByte = 0;
            curDblCount = 0;
            curTplCount = 0;
            lastInstructionId = 0;
            operCounts = OperationCounts{};

            buffer.enforceCapacityBigEnough();
            sizePosition = buffer.size;
            buffer.push64(0); // placeholder for size
            buffer.push64(B);
            buffer.size += 32; // place for operation counts
        }
        void END() { 
            buffer.push((uint8_t)Block::End);
            buffer.insert32(sizePosition, buffer.size - sizePosition);
            buffer.insert32(sizePosition + 16, operCounts.dbl);
            buffer.insert32(sizePosition + 20, operCounts.dbln);
            buffer.insert32(sizePosition + 24, operCounts.tpl);
            buffer.insert32(sizePosition + 28, operCounts.tpln);
            buffer.insert32(sizePosition + 32, operCounts.add);
            buffer.insert32(sizePosition + 36, operCounts.addn);
            buffer.insert32(sizePosition + 40, operCounts.dadd);
            buffer.insert32(sizePosition + 44, operCounts.ddbl);
        }

        // Naf
        void nafSTART(uint8_t tableSize, uint8_t initialPointIndex) {
            buffer.push((uint8_t)((uint8_t)Block::Naf | (uint8_t)(tableSize << 4)));
            if (tableSize > 1) {
                buffer.push(initialPointIndex);
                operCounts.dbl += 1;
                operCounts.add += initialPointIndex - 1;
            }
        }
        void nafSTART() { 
            nafSTART(0, 0); 
        }
        void nafEND() { 
            buffer.push((uint8_t)NafOpCode::END); 
        }
        void nafDBL() {
            buffer.push((uint8_t)NafOpCode::DBLn);
            operCounts.dbln += 1;
        }
        void nafADD(uint8_t arg) { 
            buffer.data[buffer.size-1] |= (uint8_t)NafOpCode::FullMask; 
            buffer.push((uint8_t)((uint8_t)NafOpCode::ADDn | (uint8_t)(arg << 4))); 
            // TODO: its not always correct (last op could be add/sub and not dbl)
            operCounts.addn += 1;
            operCounts.dbl += 1;
            operCounts.dbln -= 1;
        }
        void nafSUB(uint8_t arg) { 
            buffer.data[buffer.size-1] |= (uint8_t)NafOpCode::FullMask; 
            buffer.push((uint8_t)((uint8_t)NafOpCode::ADDn | (uint8_t)(arg << 4))); 
            // TODO: its not always correct (last op could be add/sub and not dbl)
            operCounts.addn += 1;
            operCounts.dbl += 1;
            operCounts.dbln -= 1;
        }
    
        // dbChain
        void dbChainSTART(uint8_t tableSize, uint8_t initialPointIndex) { 
            buffer.push((uint8_t)((uint8_t)Block::DbChain | (uint8_t)(tableSize << 4)));
            if (tableSize > 0) {
                buffer.push(initialPointIndex);
                operCounts.dbl += 1;
                operCounts.add += initialPointIndex;
            }
        }
        void dbChainSTART() {
            dbChainSTART(0, 0);
        }
        void dbChainEND() {
            if (curDblCount != 0 || curTplCount != 0) {
                dbChainAddInstruction(0, 0);
                buffer.data[lastInstructionId] |= 0x40;
                operCounts.addn -= 1;
            }
            buffer.data[lastInstructionId] |= 0x20;
            operCounts.addn -= 1;
            operCounts.add += 1;
        }
        void dbChainDBL()             { curDblCount += 1; }
        void dbChainTPL()             { curTplCount += 1; }
        void dbChainADD(uint8_t arg1) { dbChainAddInstruction(arg1, 0); }
        void dbChainSUB(uint8_t arg1) { dbChainAddInstruction(arg1, 1); }
        void dbChainAddInstruction(uint8_t arg, uint8_t signBit) {
            /*
                k - skip add/sub operation
                f - is final
                s - is subtract (instead of add)
                0kfsnnnn dddddddd
                1kfsnnnn dddddddd tttttttt
            */
            debugAssert(arg <= 15);
            uint8_t startByte = ((curTplCount > 0) << 7) | (signBit << 4) | arg;
            buffer.push(startByte);
            lastInstructionId = buffer.size - 1;
            buffer.push(curDblCount);
            if (curTplCount > 0)
                buffer.push(curTplCount);

            operCounts.addn += 1;
            if (curDblCount > 0) {
                operCounts.tpln += curTplCount;
                operCounts.dbln += curDblCount - 1;
                operCounts.dbl  += 1;
            } else {
                operCounts.tpln += curTplCount - 1;
                operCounts.tpl  += 1;
            }
            curDblCount = 0;
            curTplCount = 0;
        }


        // PRAC
        void pracSTART() { buffer.push((uint8_t)Block::Prac); }
        void pracEND()   { buffer.push((uint8_t)PracOpCode::End); }
        void pracRule(uint8_t ruleNr, bool swapBefore) {
            uint8_t byte = (uint8_t)(swapBefore << 4) | ruleNr;
            if ((pracLastByte & 0x1f) == byte && pracLastByte < 0b11100000) {
                pracLastByte = buffer.data[buffer.size - 1] += 0b00100000;
            } else {
                buffer.push(byte);
                pracLastByte = byte;
            }

            constexpr int AddCounts[] = { 3, 1, 1, 1, 1, 3, 3, 3, 1 };
            constexpr int DblCounts[] = { 0, 1, 0, 1, 1, 1, 1, 1, 1 };
            operCounts.dadd += AddCounts[ruleNr - 1];
            operCounts.ddbl += DblCounts[ruleNr - 1];
        }
    };
    void writeToFile(const std::string& fileName, Buffer& buffer) {
        std::ofstream file(fileName, std::ios::binary);
        file.write((char*)buffer.data.data(), buffer.size);
    }

    struct FileData {
        FileData() {}
        FileData(int bufferCapacity) : buffer(bufferCapacity) {}
        std::vector<uint8_t> buffer;
        OperationCounts operCounts;
    };
    std::map<uint64_t, FileData> readFromFile(const std::string& fileName) {
        std::ifstream file(fileName, std::ios::binary);
        std::map<uint64_t, FileData> result;
        uint64_t size;
        while (file.read((char*)&size, 8)) {
            FileData fileData(size - 48);
            uint64_t B;
            file.read((char*)&B, 8);
            file.read((char*)&fileData.operCounts.dbl,  4);
            file.read((char*)&fileData.operCounts.dbln, 4);
            file.read((char*)&fileData.operCounts.tpl,  4);
            file.read((char*)&fileData.operCounts.tpln, 4);
            file.read((char*)&fileData.operCounts.add,  4);
            file.read((char*)&fileData.operCounts.addn, 4);
            file.read((char*)&fileData.operCounts.dadd, 4);
            file.read((char*)&fileData.operCounts.ddbl, 4);
            file.read((char*)fileData.buffer.data(), size - 48);
            result.emplace(B, fileData);
        }
        return result;
    }
}
