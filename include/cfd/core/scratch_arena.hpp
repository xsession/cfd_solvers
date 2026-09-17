#pragma once
#include <cstddef>
#include <cstdint>
#include <new>
#include <span>
#include <stdexcept>
#include <vector>
namespace cfd::core {
class ScratchArena {
public:explicit ScratchArena(std::size_t bytes=0):storage_(bytes){}void reserve(std::size_t bytes){if(offset_!=0)throw std::logic_error("cannot reserve active scratch arena");storage_.resize(bytes);}void reset()noexcept{offset_=0;}[[nodiscard]] std::size_t capacity()const noexcept{return storage_.size();}[[nodiscard]] std::size_t used()const noexcept{return offset_;}
    template<class T> std::span<T> allocate(std::size_t count){const std::size_t align=alignof(T);std::size_t pos=(offset_+align-1)&~(align-1);const std::size_t bytes=count*sizeof(T);if(pos+bytes>storage_.size())throw std::bad_alloc();auto* ptr=reinterpret_cast<T*>(storage_.data()+pos);offset_=pos+bytes;for(std::size_t i=0;i<count;++i)::new(static_cast<void*>(ptr+i))T{};return {ptr,count};}
private:std::vector<std::byte>storage_;std::size_t offset_{};
};
} // namespace cfd::core
