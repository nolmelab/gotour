package test

import (
	"testing"
	"unsafe"
)

// https://go101.org/article/unsafe.html
// Learn about unsafe pointers
//

func TestIntro(t *testing.T) {
	// type-safe pointers vs. type-unsafe pointers
	// unsafe.Pointer is like void* in C

	var x struct {
		a int64
		b bool
		c string
	}

	const M, N = unsafe.Sizeof(x.c), unsafe.Sizeof(x)
	t.Log(M, N)

	t.Log(unsafe.Alignof(x.a))
	t.Log(unsafe.Alignof(x.b))
	t.Log(unsafe.Alignof(x.c))

	t.Log(unsafe.Offsetof(x.a))
	t.Log(unsafe.Offsetof(x.b))
	t.Log(unsafe.Offsetof(x.b))
}

func TestIntroEmbed(t *testing.T) {
	type T struct {
		c string
	}

	type S struct {
		b bool
	}

	var x struct {
		a int64
		*S
		T
	}

	t.Log(unsafe.Offsetof(x.a))
	t.Log(unsafe.Offsetof(x.S))
	t.Log(unsafe.Offsetof(x.T))

	t.Log(unsafe.Offsetof(x.c))
	// Pointer S is embedded with a pointer
	// t.Log(unsafe.Offsetof(x.b))
	t.Log(unsafe.Offsetof(x.S.b))
}

func TestIntroPointerArith(t *testing.T) {
	a := [16]int{3: 3, 9: 9, 11: 11}
	t.Log(a)
	eleSize := int(unsafe.Sizeof(a[0]))
	t.Log(eleSize)
	p9 := &a[9]
	up9 := unsafe.Pointer(p9)
	p3 := (*int)(unsafe.Add(up9, -6*eleSize))
	t.Log(*p3)

	s := unsafe.Slice(p9, 5)[:3]
	t.Log(s)
	t.Log(len(s), cap(s))

	u := unsafe.Slice((*int)(nil), 0)
	t.Log(u == nil)
}

func TestIntroByteSlice(t *testing.T) {
	s := "Hello"
	bs := unsafe.Slice(unsafe.StringData(s), len(s))
	s2 := unsafe.String(unsafe.SliceData(bs), len(bs))

	t.Log(s == s2)
	t.Log(s2)
}

func TestPointerConversion(t *testing.T) {
	var v int16 = 7200
	p1 := unsafe.Pointer(&v)
	p2 := unsafe.Add(p1, 1)
	bp := (*int8)(p2)
	v2 := *bp

	t.Log(v2)
	t.Log(v >> 8)

	// 함수를 사용한다는 점 외에는 C와 거의 같다.
}

func TestFact1(t *testing.T) {
	// Fact 1. unsafe pointers are pointers and uintptr values are integers
	// var i uintptr // builtin
	// i = 3
	// p1 := unsafe.Pointer(i)

	// segmentation fault as expected
	// v := (*int)(p1)
	// t.Log(*v)
}

func TestFact2(t *testing.T) {
	// TestFact1에서 이미 uintptr의 위험성을 보였다.
	// uintptr에서 unsafe.Pointer로 변환이 가능하다.
	// uintptr은 int처럼 연산이 가능하다.
}

func TestFact3(t *testing.T) {
	// Fact 3. the addresses of some values might change at run time
	// 스택의 크기 변화나 다른 이유로 인해 메모리 주소가 바뀔 수 있다.
}

func TestFact4(t *testing.T) {
	// Fact 4. the life range of a value at run time may be not as large as it looks in code
	type T struct {
		x int
		y *[1 << 23]byte
	}

	// v := T{y: new([1<<23]byte)}
	// p := uintptr(unsafe.Pointer(&v.y[0]))

	// uintptr은 GC에서 고려하지 않는다.
	// 따라서, v.y 접근은 위험할 수 있다. 진짜? 컴파일러 최적화로 가능하다.
}

func TestFact5(t *testing.T) {
	// Fact 5. *unsafe.Pointer is a general safe pointer type
	x := 123
	p := unsafe.Pointer(&x)
	pp := &p                  // *unsafe.Pointer
	p = unsafe.Pointer(pp)    // p is again of unsafe.Pointer type
	pp = (*unsafe.Pointer)(p) // pp is again of *unsafe.Pointer type

	t.Log(*pp) // pp holds the pointer(address) to p
}

func TestSafePattern1(t *testing.T) {
	// Pattern 1. convert a *T1 to unsafe Pointer, then convert the unsafe pointer
	// value to *T2

	type MyString string
	ms := []MyString{"C", "C++", "Go"}
	t.Log(ms)

	ss := *(*[]string)(unsafe.Pointer(&ms))
	ss[1] = "Zig"
	t.Log(ss)

	ms = *(*[]MyString)(unsafe.Pointer(&ss))
	t.Log(ms)

	// Slice can be used also

	// 변환은 주의해서 사용해야 한다.
}

func TestSafePattern2(t *testing.T) {
	// Pattern 2. convert unsafe pointer to uintptr, then use the uinptr value

	type T struct {
		a int
	}

	var v T
	t.Logf("%p\n", &v)
	t.Logf("0x%x\n", uintptr(unsafe.Pointer(&v)))
}

func TestSafePattern3(t *testing.T) {
	type T struct {
		x bool
		y [3]int16
	}

	const N = unsafe.Offsetof(T{}.y)
	const M = unsafe.Sizeof(T{}.y[0])

	v := T{y: [3]int16{123, 456, 789}}
	p := unsafe.Pointer(&v)
	ty2 := (*int16)(unsafe.Pointer(uintptr(p) + N + M + M))
	t.Log(*ty2)

	// comment: GC나 메모리 변경으로 인한 주소 변경 때문에 주로 힙에서 할당하고
	// 누군가 들고 있는 포인터에 대해서만 unsafe를 사용하는 것이 좋다.
	// 스택 변수도 GC 될 수 있고 주소가 바뀔 수도 있다.
}

func TestSafePattern4(t *testing.T) {

}

func TestUnsafeSwapBytes(t *testing.T) {

}
