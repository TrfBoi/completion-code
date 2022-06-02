// Package main 程序主包
package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io/fs"
	"io/ioutil"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
)

func main() {
	// 1. 制作基因表达输入矩阵
	//s genMatrix()
	// 2. 制作生存分析输入矩阵
	// genSurvivalInputMatrix()
}

// 压缩counts文件名切片
var gzFiles []string
var unGZFilePath = "D:/Users/Lenovo/Desktop/毕设/rnaseq_data/"

// func init() {
// 	err := os.Mkdir("D:/Users/Lenovo/Desktop/毕设/rnaseq_data", os.ModePerm)
// 	if err != nil {
// 		panic("创建解压目录失败")
// 	}
// }

// genMatrix 生成基因-样本表达量矩阵
func genMatrix() {
	// 1. 解压并集中文件
	// countsPath := "D:/Users/Lenovo/Desktop/毕设/gdc_download_20220227_071028.689186"
	// err := searchAllGZFileFromDir(countsPath)
	// if err != nil {
	// 	panic("RNA-Seq数据目录输入出错")
	// }
	// for _, gzFileName := range gzFiles {
	// 	gzFileFD, err := os.Open(gzFileName)
	// 	if err != nil {
	// 		panic("GZ格式文件打开出错" + gzFileName)
	// 	}
	// 	defer gzFileFD.Close()
	// 	gzReader, err := gzip.NewReader(gzFileFD)
	// 	if err != nil {
	// 		panic("读取GZ文件内容出错" + gzFileName)
	// 	}
	// 	defer gzReader.Close()
	// 	unGZFile, err := os.Create(unGZFilePath + gzReader.Header.Name)
	// 	if err != nil {
	// 		panic("创建文件出错")
	// 	}
	// 	defer unGZFile.Close()
	// 	buf, err := ioutil.ReadAll(gzReader)
	// 	if err != nil {
	// 		panic("读取压缩文件失败")
	// 	}
	// 	unGZFile.Write(buf)
	// }
	// 2. 生成表达量矩阵geneSampleMatrix
	geneSampleMatrix, fileIds := doGenMatrix()
	// 2.1 利用metadata中文件名与样本名对应关系修改列名
	filesId2SampleId(fileIds)
	// 2.2 geneId转换为为基因名
	// geneId2GeneName(geneSampleMatrix)
	// 3 输出矩阵到文件
	os.Remove("D:/Users/Lenovo/Desktop/毕设/matrix.txt")
	matrixFd, err := os.Create("D:/Users/Lenovo/Desktop/毕设/matrix.txt")
	if err != nil {
		panic(err)
	}
	defer matrixFd.Close()
	matrixFd.WriteString("geneName")
	for _, geneId := range fileIds {
		// if i == 0 {
		// 	matrixFd.WriteString(geneId)
		// 	continue
		// }
		matrixFd.WriteString(" " + geneId)
	}
	matrixFd.WriteString("\n")
	for _, gene := range geneSampleMatrix {
		for geneId, counts := range gene {
			// 去掉基因ID版本号
			geneId = strings.Split(geneId, ".")[0]
			matrixFd.WriteString(geneId + " ")
			matrixFd.WriteString(counts + "\n")
		}
	}
}

// searchAllGZFileFromDir 找出指定目录下的所有gz格式压缩文件名
func searchAllGZFileFromDir(dirName string) error {
	err := filepath.Walk(dirName, func(path string, info fs.FileInfo, err error) error {
		if info.IsDir() {
			return nil
		}
		if len(path) <= 1 || path[len(path)-2:] != "gz" {
			return nil
		}
		gzFiles = append(gzFiles, path)
		return nil
	})
	fmt.Printf("样本数: %d \n", len(gzFiles))
	return err
}

// doGenMatrix 合并解压后的count文件
func doGenMatrix() ([]map[string]string, []string) {
	geneSampleMatrix := make([]map[string]string, 60488)
	fileIds := make([]string, 151)
	idIndex := 0
	filepath.Walk(unGZFilePath, func(path string, info fs.FileInfo, err error) error {
		if info.IsDir() {
			return nil
		}
		fd, err := os.Open(path)
		if err != nil {
			panic("打开解压缩后的count文件失败")
		}
		defer fd.Close()
		fileIds[idIndex] = info.Name() + ".gz"
		readLineFd := bufio.NewScanner(fd)
		i := 0
		for {
			if !readLineFd.Scan() {
				break
			}
			line := readLineFd.Text()
			curGeneCount := strings.Split(line, "	")
			if len(geneSampleMatrix[i]) == 0 {
				geneSampleMatrix[i] = map[string]string{curGeneCount[0]: curGeneCount[1]}
				i++
				continue
			}
			geneSampleMatrix[i][curGeneCount[0]] = geneSampleMatrix[i][curGeneCount[0]] + " " + curGeneCount[1]
			i++
		}
		idIndex++
		return nil
	})
	return geneSampleMatrix, fileIds
}

// filesId2SampleId 将文件id转换为文件对应的样本id
func filesId2SampleId(fileIds []string) {
	metadataFilePath := "D:/Users/Lenovo/Desktop/毕设/metadata.cart.2022-02-27.json"
	fd, err := os.Open(metadataFilePath)
	if err != nil {
		panic(err)
	}
	defer fd.Close()
	jsonByte, err := ioutil.ReadAll(fd)
	if err != nil {
		panic(err)
	}
	metadata := make([]map[string]interface{}, 1)
	err = json.Unmarshal(jsonByte, &metadata)
	if err != nil {
		panic(err)
	}
	for i, fileId := range fileIds {
		for _, meta := range metadata {
			metaFileId := meta["file_name"].(string)
			if metaFileId == fileId {
				associatedEntities := meta["associated_entities"].([]interface{})
				fileIds[i] = associatedEntities[0].(map[string]interface{})["entity_submitter_id"].(string)
				break
			}
		}
	}
}

// geneId2GeneName 基因id转换为对应基因名, 方便人
func geneId2GeneName(geneSampleMatrix []map[string]string) {
	// 1. 打开gtf格式基因注释文件
	fd, err := os.Open("D:/Users/Lenovo/Desktop/毕设/Homo_sapiens.GRCh38.105.chr.gtf")
	if err != nil {
		panic(err)
	}
	defer fd.Close()
	// 2. 双重遍历依次替换基因id
	done := make([]bool, 60488)
	doneCount := 0
	matchGeneIdName := regexp.MustCompile("gene_id \"(.+?)\"; gene_version \"(.+?)\"; gene_name \"(.+?)\";")
	scan := bufio.NewScanner(fd)
	for {
		if !scan.Scan() {
			break
		}
		line := scan.Text()
		// 不是gene类型行, 舍去
		if !matchGeneIdName.MatchString(line) {
			continue
		}
		for i, gene := range geneSampleMatrix {
			// 转换过的快速进入下一个
			if done[i] {
				continue
			}
			for key, value := range gene {
				geneInfos := matchGeneIdName.FindStringSubmatch(line)
				newKey := strings.Split(key, ".")[0]
				if geneInfos[1] == newKey {
					geneSampleMatrix[i][geneInfos[3]] = value
					delete(geneSampleMatrix[i], key)
					done[i] = true
					doneCount++
				}
			}
			// 一行只会匹配一个(待验证)
			if done[i] {
				break
			}
		}
	}
	fmt.Println(geneSampleMatrix)
	fmt.Println(doneCount)
	fmt.Println(len(geneSampleMatrix))
}

// genSurvivalInputMatrix 生成生存分析输入矩阵
func genSurvivalInputMatrix() {
	// 1.打开临床json格式文件并解析
	fd, err := os.Open("D:/Users/Lenovo/Desktop/毕设/clinical.cart.2022-02-27.json")
	if err != nil {
		panic(err)
	}
	defer fd.Close()
	clinicalByte, err := ioutil.ReadAll(fd)
	if err != nil {
		panic(err)
	}
	clinical := make([]map[string]interface{}, 151)
	err = json.Unmarshal(clinicalByte, &clinical)
	if err != nil {
		panic(err)
	}
	// 2.取出生存时间, 是否死亡, 性别, 年龄等信息
	survivalInputMatrix := make(map[string][]string, 5)
	for i, patient := range clinical {
		diagnoses := patient["diagnoses"].([]interface{})[0].(map[string]interface{})

		id := strings.Split(diagnoses["submitter_id"].(string), "_")[0]
		if len(survivalInputMatrix["id"]) == 0 {
			survivalInputMatrix["id"] = []string{id}
		} else {
			survivalInputMatrix["id"] = append(survivalInputMatrix["id"], id)
		}
		demographic := patient["demographic"].(map[string]interface{})
		if len(survivalInputMatrix["alive"]) == 0 {
			survivalInputMatrix["alive"] = []string{"1"}
		} else {
			survivalInputMatrix["alive"] = append(survivalInputMatrix["alive"], "1")
		}
		gender := demographic["gender"].(string)
		if len(survivalInputMatrix["gender"]) == 0 {
			survivalInputMatrix["gender"] = []string{gender}
		} else {
			survivalInputMatrix["gender"] = append(survivalInputMatrix["gender"], gender)
		}
		race := demographic["race"].(string)
		if len(survivalInputMatrix["race"]) == 0 {
			survivalInputMatrix["race"] = []string{race}
		} else {
			survivalInputMatrix["race"] = append(survivalInputMatrix["race"], race)
		}
		tmpAge := demographic["age_at_index"].(float64)
		age := strconv.FormatFloat(tmpAge, 'f', 1, 64)
		if len(survivalInputMatrix["age_at_index"]) == 0 {
			survivalInputMatrix["age_at_index"] = []string{age}
		} else {
			survivalInputMatrix["age_at_index"] = append(survivalInputMatrix["age_at_index"], age)
		}
		var time float64
		ok := true
		if demographic["vital_status"].(string) == "Alive" {
			time, ok = diagnoses["days_to_last_follow_up"].(float64)
		} else {
			survivalInputMatrix["alive"][i] = "0"
			time, ok = demographic["days_to_death"].(float64)
		}
		if !ok {
			time = -1
		}
		timeString := strconv.FormatFloat(float64(time), 'f', 1, 64)
		if len(survivalInputMatrix["time"]) == 0 {
			survivalInputMatrix["time"] = []string{timeString}
		} else {
			survivalInputMatrix["time"] = append(survivalInputMatrix["time"], timeString)
		}
	}
	// 3. 校验并存入文件
	os.Remove("D:/Users/Lenovo/Desktop/毕设/survivalInputMatrix.txt")
	out, err := os.Create("D:/Users/Lenovo/Desktop/毕设/survivalInputMatrix.txt")
	if err != nil {
		panic(err)
	}
	defer out.Close()
	for key, value := range survivalInputMatrix {
		if len(value) != 151 {
			panic("生成survivalInputMatrix失败, 数据数目不对")
		}
		out.WriteString(key)
		for _, s := range value {
			out.WriteString("\t" + s)
		}
		out.WriteString("\n")
	}
}
