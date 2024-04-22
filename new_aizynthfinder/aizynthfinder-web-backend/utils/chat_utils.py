import os
import sys
import torch
from peft import PeftModel
from rdkit import Chem
from rdkit.Chem import MolFromSmiles
from selfies import decoder
from transformers import GenerationConfig, LlamaForCausalLM, LlamaTokenizer
from utils.prompter import Prompter
from translate import Translator
translator = Translator(to_lang="en")
if torch.cuda.is_available():
    device = "cuda"
else:
    device = "cpu"

try:
    if torch.backends.mps.is_available():
        device = "mps"
except:  # noqa: E722
    pass


def main(
        CLI: bool = False,
        protein: bool = False,
        load_8bit: bool = False,
        base_model: str = "",
        lora_weights: str = "zjunlp/llama-molinst-molecule-7b",
        prompt_template: str = "",
        server_name: str = "0.0.0.0",  # Allows to listen on all interfaces by providing '0.
        share_gradio: bool = False,
):
    base_model = base_model or os.environ.get("BASE_MODEL", "")
    assert (
        base_model
    ), "Please specify a --base_model, e.g. --base_model='decapoda-research/llama-7b-hf'"

    prompter = Prompter(prompt_template)
    if protein == False:
        tokenizer = LlamaTokenizer.from_pretrained(base_model)
    else:
        tokenizer = LlamaTokenizer.from_pretrained(base_model, bos_token='<s>', eos_token='</s>', add_bos_token=True,
                                                   add_eos_token=False)
    if device == "cuda":
        model = LlamaForCausalLM.from_pretrained(
            base_model,
            load_in_8bit=load_8bit,
            torch_dtype=torch.float16,
            # device_map="auto",
            device_map={"": 0}
        )
        if protein == False:
            model = PeftModel.from_pretrained(
                model,
                lora_weights,
                torch_dtype=torch.float16,
                device_map={"": 0},
            )
    elif device == "mps":
        model = LlamaForCausalLM.from_pretrained(
            base_model,
            device_map={"": device},
            torch_dtype=torch.float16,
        )
        if protein == False:
            model = PeftModel.from_pretrained(
                model,
                lora_weights,
                device_map={"": device},
                torch_dtype=torch.float16,
            )
    else:
        model = LlamaForCausalLM.from_pretrained(
            base_model, device_map={"": device}, low_cpu_mem_usage=True
        )
        if protein == False:
            model = PeftModel.from_pretrained(
                model,
                lora_weights,
                device_map={"": device},
            )

    # unwind broken decapoda-research config
    model.config.pad_token_id = tokenizer.pad_token_id = 0  # unk
    model.config.bos_token_id = 1
    model.config.eos_token_id = 2

    if not load_8bit:
        model.half()  # seems to fix bugs for some users.

    model.eval()
    if torch.__version__ >= "2" and sys.platform != "win32":
        model = torch.compile(model)

    def evaluate(
            instruction,
            input=None,
            temperature=0.1,
            top_p=0.75,
            top_k=40,
            num_beams=4,
            repetition_penalty=1,
            max_new_tokens=128,
            **kwargs,
    ):

        prompt = prompter.generate_prompt(instruction, input)
        inputs = tokenizer(prompt, return_tensors="pt")
        input_ids = inputs["input_ids"].to(device)
        if protein == False:
            do_sample = False
        else:
            do_sample = True
        generation_config = GenerationConfig(
            do_sample=do_sample,
            temperature=temperature,
            top_p=top_p,
            top_k=top_k,
            num_beams=num_beams,
            repetition_penalty=repetition_penalty,
            **kwargs,
        )
        with torch.no_grad():
            generation_output = model.generate(
                input_ids=input_ids,
                generation_config=generation_config,
                return_dict_in_generate=True,
                output_scores=True,
                max_new_tokens=max_new_tokens,
            )
        s = generation_output.sequences[0]
        output = tokenizer.decode(s)
        re = prompter.get_response(output)
        # remove the last ']' or ‘#’
        last_bracket_index = re.find('#')
        if last_bracket_index != -1:
            re = re[:last_bracket_index]

        last_bracket_index = re.rfind(']')
        if last_bracket_index != -1:
            re = re[:last_bracket_index + 1]

        return re

    # 定义 evaluate 函数
    def evaluate_instruction(instruction, input_text, temperature=0.1, top_p=0.75, top_k=40, num_beams=4,
                             repetition_penalty=1, max_new_tokens=128):
        # 如果 instruction 是中文，翻译成英文
        if is_chinese(instruction):
            instruction = translator.translate(instruction)
        output = evaluate(instruction, input=input_text, temperature=temperature, top_p=top_p, top_k=top_k,
                          num_beams=num_beams, repetition_penalty=repetition_penalty, max_new_tokens=max_new_tokens)
        # 如果 output 是英文，翻译成中文
        if not is_chinese(output):
            output = Translator(to_lang="zh").translate(output)
        return output

    def is_chinese(text):
        for char in text:
            if '\u4e00' <= char <= '\u9fff':
                return True
        return False

    def is_smiles(text):
        try:
            MolFromSmiles(text)
            return True
        except Exception as e:
            return False

    def is_selfies(text):
        try:
            decoded = decoder(text)
            return True
        except Exception as e:
            return False






