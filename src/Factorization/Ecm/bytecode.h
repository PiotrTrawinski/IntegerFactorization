#pragma once
#include <cstdint>
#include <vector>

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

    struct Buffer {
        std::vector<uint8_t> data;
        int size = 0;

        Buffer() {
            data.resize(100'000);
        }
        void push(uint8_t n) {
            data[size++] = n;
        }
        void enforceCapacityBigEnough() {
            if (size + 10'000 >= data.capacity()) {
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

        Reader(const std::vector<uint8_t>& buffer) : buffer(buffer) {}

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

    struct Writer {
        Buffer buffer;
        uint8_t pracLastByte = 0;
        int dblCount = 0;
        int tplCount = 0;
        int lastInstructionId = 0;

        // Naf
        void END()                 { buffer.push((uint8_t)Block::End); }
        void nafSTART(uint8_t tableSize, uint8_t initialPointIndex) {
            buffer.push((uint8_t)((uint8_t)Block::Naf | (uint8_t)(tableSize << 4)));
            if (tableSize > 1) {
                buffer.push(initialPointIndex);
            }
        }
        void nafSTART()            { nafSTART(0, 0); }
        void nafEND()              { buffer.push((uint8_t)NafOpCode::END); }
        void nafDBL()              { buffer.push((uint8_t)NafOpCode::DBLn); }
        void nafADD(uint8_t arg)   { buffer.data[buffer.size-1] |= (uint8_t)NafOpCode::FullMask; buffer.push((uint8_t)((uint8_t)NafOpCode::ADDn | (uint8_t)(arg << 4))); }
        void nafSUB(uint8_t arg)   { buffer.data[buffer.size-1] |= (uint8_t)NafOpCode::FullMask; buffer.push((uint8_t)((uint8_t)NafOpCode::ADDn | (uint8_t)(arg << 4))); }
    
        // dbChain
        void dbChainSTART(uint8_t tableSize, uint8_t initialPointIndex) { 
            buffer.push((uint8_t)((uint8_t)Block::DbChain | (uint8_t)(tableSize << 4)));
            if (tableSize > 1) {
                buffer.push(initialPointIndex);
            }
        }
        void dbChainSTART() {
            dbChainSTART(0, 0);
        }
        void dbChainEND() {
            if (dblCount != 0 || tplCount != 0) {
                dbChainAddInstruction(0, 0);
                buffer.data[lastInstructionId] |= 0x40;
            }
            buffer.data[lastInstructionId] |= 0x20;
        }
        void dbChainDBL()             { dblCount += 1; }
        void dbChainTPL()             { tplCount += 1; }
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
            uint8_t startByte = ((tplCount > 0) << 7) | (signBit << 4) | arg;
            buffer.push(startByte);
            lastInstructionId = buffer.size - 1;
            buffer.push(dblCount);
            if (tplCount > 0)
                buffer.push(tplCount);

            dblCount = 0;
            tplCount = 0;
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
        }
    };

    void writeToFile(std::ofstream& file, Buffer& buffer) {
        int32_t size = buffer.size;
        file.write((char*)&size, 4);
        file.write((char*)buffer.data.data(), buffer.size);
    }
}
