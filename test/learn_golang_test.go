package test

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestMultipleAssignment(t *testing.T) {
	i := 10
	j := 20

	i, j = i+1, i // 값으로 전달되는 것과 같다. 즉, i=i+1이 반영된 i 값이 아니라 이전 값을 사용
	t.Log(i)
	t.Log(j)

	assert.True(t, i == 11 && j == 10)
}
